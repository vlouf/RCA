# Python standard library
import os
import traceback

# Other libraries
import pyart
import cftime
import netCDF4
import numpy as np


def _read_with_pyart(infile, dbz_name, zdr_name, rhohv_name):
    """
    Parameters:
    ===========
    infile: str
        Radar file name.
    dbz_name: str
        Radar reflectivity field name.
    zdr_name: str
        Radar differential reflectivity field name.

    Returns:
    ========
    volume_date: datetime
        Datetime for input volume.
    r: ndarray
        Radar range.
    azi: ndarray
        Radar azimuth.
    reflec: ndarray
        Radar reflectivity.
    zdr: ndarray
        Radar differential reflectivity.
    """
    # Check file extension.
    file_extension = os.path.splitext(infile)[-1]
    try:
        if file_extension == ".h5" or file_extension == ".H5":
            radar = pyart.aux_io.read_odim_h5(infile, file_field_names=True)
        else:
            radar = pyart.io.read(infile)
    except Exception:
        print("Could not read input file", os.path.basename(infile))
        return None

    volume_date = cftime.num2date(radar.time['data'][0],
                                  radar.time['units'],
                                  only_use_cftime_datetimes=False,
                                  only_use_python_datetimes=True)

    # Extract first elevation only
    rslice = radar.get_slice(0)
    # Extract range/azimuth
    r = radar.range['data'].astype(int)
    azi = radar.azimuth['data'][rslice]

    # Get reflectivity
    reflec = None
    for dname in [dbz_name, 'reflectivity', 'DBZ', 'DBZH']:
        try:
            reflec = radar.fields[dname]['data'][rslice].filled(np.NaN)
            break
        except KeyError:
            continue

    if reflec is None:
        print("Wrong RHOHV/DBZ field names provided. The field names in radar files are:")
        raise KeyError("Wrong field name provided")

    try:
        if zdr_name is not None:
            zdr = radar.fields[zdr_name]['data'][rslice].filled(np.NaN)
        else:
            zdr = None
    except KeyError:
        print("Wrong RHOHV/DBZ field names provided. The field names in radar files are:")
        raise KeyError("Wrong field name provided")

    try:
        rhohv = radar.fields[rhohv_name]['data'][rslice]
    except Exception:
        # print("Problem with cross-correlation ratio field. Maybe missing? Continuing without it.")
        rhohv = None

    try:  # In case rhohv is a MaskedArray
        rhohv = rhohv.filled(np.NaN)
    except AttributeError:
        pass

    return volume_date, r, azi, reflec, zdr, rhohv


def _read_with_netcdf(infile, dbz_name, zdr_name, rhohv_name):
    with netCDF4.Dataset(infile, "r") as ncid:
        # Extract datetime
        volume_date = cftime.num2date(ncid['time'][0], ncid['time'].units,
                                      only_use_cftime_datetimes=False,
                                      only_use_python_datetimes=True)
        # Get first sweep
        sweep = ncid["sweep_start_ray_index"][:]
        stsw = sweep[0]
        edsw = sweep[1] - 1

        # Extract range and azimuth
        azi = ncid["azimuth"][stsw:edsw]
        r = ncid["range"][:]

        # Get reflectivity and ZDR.
        try:
            refl = ncid[dbz_name][stsw:edsw, :].filled(np.NaN)
            if zdr_name is not None:
                zdr = ncid[zdr_name][stsw:edsw, :].filled(np.NaN)
            else:
                zdr = None
        except KeyError:
            print("Wrong RHOHV/DBZ field names provided. The field names in radar files are:")
            print(radar.fields.keys())
            raise KeyError("Wrong field name provided")

        # Extract RHOHV
        try:
            rhohv = ncid[rhohv_name][stsw:edsw, :]
        except Exception:
            traceback.print_exc()
            print("Problem with cross-correlation ratio field. Maybe missing? Continuing without it.")
            rhohv = None

        try:  # In case rhohv is a MaskedArray
            rhohv = rhohv.filled(np.NaN)
        except AttributeError:
            pass

    return volume_date, r, azi, refl, zdr, rhohv


def read_data(infile, dbz_name="DBZ", zdr_name=None, rhohv_name="RHOHV"):
    """
    Parameters:
    ===========
    infile: str
        Radar file name.
    dbz_name: str
        Radar reflectivity field name.
    zdr_name: str
        Radar differential reflectivity field name.

    Returns:
    ========
    volume_date: datetime
        Datetime for input volume.
    r: ndarray
        Radar range.
    azi: ndarray
        Radar azimuth.
    reflec: ndarray
        Radar reflectivity.
    zdr: ndarray
        Radar differential reflectivity.
    """
    # Check file extension.
    file_extension = os.path.splitext(infile)[-1]

    if file_extension == ".nc" or file_extension == ".NC":
        rslt = _read_with_netcdf(infile, dbz_name, zdr_name, rhohv_name)
    else:
        rslt = _read_with_pyart(infile, dbz_name, zdr_name, rhohv_name)

    if rslt is None:
        volume_date, r, azi, reflec, zdr, rhohv = [None] * 6
    else:
        volume_date, r, azi, reflec, zdr, rhohv = rslt

    return volume_date, r, azi, reflec, zdr, rhohv


def write_ncfile(outfilename, xdate, rca, rca_zdr, rain, gnrl_meta):
    """
    Write data to netCDF4 file.

    Parameters:
    ===========
    outfilename: str
        Output file name saving mask.
    xdate: np.array
        Time dimension.
    rca: np.array
        Ground clutter 95 percentile reflectivity value.
    gnrl_meta: dict
        Metadata dictionnary
    """
    dim_len = len(xdate)
    time_units = "seconds since 1970-01-01"
    date = xdate.tolist()
    time = [netCDF4.date2num(d, time_units) for d in date]

    # Write netCDF4 file.
    with netCDF4.Dataset(outfilename, "w", format="NETCDF4") as rootgrp:
        # Create dimension
        rootgrp.createDimension("time", dim_len)

        # Create variables.
        ncr = rootgrp.createVariable('time', 'f8', ("time",), zlib=True)
        nca = rootgrp.createVariable('rca', 'f8', ("time",), zlib=True)
        ncrain = rootgrp.createVariable('rain', 'f8', ("time",), zlib=True)

        # Assign values.
        ncr[:] = time
        nca[:] = rca
        ncrain[:] = rain

        # Set units.
        ncr.units = time_units
        nca.units = "dBZ"
        ncrain.units = "mm/h"

        if len(rca_zdr) > 0:
            ncz = rootgrp.createVariable('rca_zdr', 'f8', ("time",), zlib=True)
            ncz[:] = rca_zdr
            ncz.units = "dB"

        # Set main metadata
        for mykey in gnrl_meta.keys():
            rootgrp.setncattr_string(mykey, gnrl_meta[mykey])

    return None
