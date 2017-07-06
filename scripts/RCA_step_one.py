"""
First part of the ground-clutter-monitoring calibration technique. This part
handles the creation of the clutter mask.

@title: RCA_step_one
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Bureau of Meteorology/Monash University
@date: 5/07/2017
@version: 0.99

.. autosummary::
    :toctree: generated/

    plot_freq_map
    write_ncfile
    multproc_buffer_create_clut_map
    main
"""
# Python standard library
import os
import hashlib
import datetime
import warnings
import traceback

from multiprocessing import Pool

# Other modules. Matplotlib must be imported first
import matplotlib
matplotlib.use('Agg')  # <- Reason why matplotlib is imported first.
import matplotlib.pyplot as pl
import pyart
import netCDF4
import numpy as np

# Custom modules.
from processing_codes import cmask_code
from processing_codes.raijin_tools import get_files


def plot_freq_map(outfilename_fig, rrange, freq):
    """
    Plot the clutter frequency map.

    Parameters:
    ===========
        outfilename_fig: str
            Output file name for figure.
        rrange: np.array
            Radar range
        freq: np.array
            Clutter frequency array.
    """
    def azi2angl(azi):
        """
        Convert azimuth angle to true angle (in degree).
        """
        h = 450 - azi
        h[h > 360] = h[h > 360] - 360
        return h

    # Mask useless data.
    freq = np.ma.masked_where(freq == 0, freq)

    # Create azimuth dimension
    r = rrange
    aziclut = np.linspace(0, 360, freq.shape[0], endpoint=False)
    theta = azi2angl(aziclut)

    # Get cartesians dimension.
    [TH, R] = np.meshgrid(theta, r, indexing='ij')
    x = R*np.cos(TH*np.pi/180)
    y = R*np.sin(TH*np.pi/180)

    pl.figure()
    pl.pcolor(x, y, freq, cmap='pyart_Theodore16', vmin=0, vmax=100)
    for r in [5e3, 10e3]:
        pl.plot(r*np.cos(theta*np.pi/180), r*np.sin(theta*np.pi/180), 'k-')
    pl.axis('square')
    pl.axis((-10e3, 10e3, -10e3, 10e3))
    pl.title("Non-meteorological echoes")
    pl.xlabel("Distance from radar origin (m)")
    pl.ylabel("Distance from radar origin (m)")
    c0 = pl.colorbar()
    c0.set_label("Frequency (%)")
    pl.savefig(outfilename_fig, dpi=250)
    pl.close()

    return None


def write_ncfile(outfilename, clutter_r, clutter_azi, gnrl_meta):
    """
    Write data to netCDF4 file.

    Parameters:
    ===========
        outfilename: str
            Output file name saving mask.
        clutter_r: np.array
            Array of clutter range position.
        clutter_azi: np.array
            Array of clutter azimuthal position.
        gnrl_meta: dict
            Metadata dictionnary
    """
    # Check that range and azimuth are the same length.
    if len(clutter_r) != len(clutter_azi):
        print("Invalid dimension.")
        return None

    dim_len = len(clutter_r)

    # Write netCDF4 file.
    with netCDF4.Dataset(outfilename, "w", format="NETCDF4") as rootgrp:
        # Create dimension
        rootgrp.createDimension("x", dim_len)

        # Create variables.
        ncr = rootgrp.createVariable('range', 'i8', ("x",), zlib=True)
        nca = rootgrp.createVariable('azimuth', 'f8', ("x",), zlib=True)

        # Assign values.
        ncr[:] = clutter_r
        nca[:] = clutter_azi

        # Set units.
        ncr.units = "meters"
        nca.units = "degrees"

        # Set main metadata
        for mykey in gnrl_meta.keys():
            rootgrp.setncattr_string(mykey, gnrl_meta[mykey])

    return None


def multproc_buffer_create_clut_map(infile):
    """
    Buffer function for multiprocessing and handleing errors. No actual
    computation is done here. It just calls the function that will extract the
    clutter range and azimuth position.

    Parameters:
    ===========
        infile: str
            Radar file name.

    Returns:
    ========
        r_clutt: np.array
            Range position of non-meteorological echoes.
        azi_clutt: np.array
            Azimuth position of non-meteorological echoes.
    """
    try:
        radar = pyart.io.read(infile)
    except Exception:
        print("Could not read input file", os.path.basename(infile))
        return None

    try:
        r_clutt, azi_clutt = cmask_code.get_clutter_position(radar,
                                                             dbz_name=DBZ_FIELD_NAME,
                                                             rhohv_name=RHOHV_FIELD_NAME,
                                                             refl_thrld=45,
                                                             rhohv_thrld=0.6,
                                                             maxrange=10e3)
    except Exception:
        print("Problem with this file:", os.path.basename(infile))
        traceback.print_exc()
        return None

    return r_clutt, azi_clutt


def main():
    flist = get_files(INPUT_DIR)
    # Create the name of output files (figure and mask).
    outfilename_suffix = hashlib.md5(" ".join(flist).encode('utf-8')).hexdigest()[:16]
    # netCDF4 File
    outfilename_save = "CLUTTER_map_{}.nc".format(outfilename_suffix)
    outfilename_save = os.path.join(OUTPUT_DIR, outfilename_save)
    # PNG figure.
    outfilename_fig = "Frequency_map_{}.png".format(outfilename_suffix)
    outfilename_fig = os.path.join(OUTPUT_DIR, outfilename_fig)

    if os.path.isfile(outfilename_save):
        print("Output data file already exists. Doing nothing.")
        return None

    # Extract info of range and azimuth from one file
    radar = pyart.io.read(flist[0])
    rrange = radar.range['data'].astype(int)
    azimuth = radar.azimuth['data'][radar.get_slice(0)]
    nbfile = len(flist)

    if nbfile < 100:
        print("Only {} file(s) found.".format(nbfile))
        print("Need more data to generate a valid clutter mask. Doing nothing.")
        return None

    # Start multiprocessing
    with Pool(NCPU) as pool:
        rslt = pool.map(multproc_buffer_create_clut_map, flist)

    print("Non-meteorological echoes extracted.")
    # Unpack multiprocessing rslt.
    range_tot = np.array([], dtype=int)
    azi_tot = np.array([])
    for rslice, azislice in rslt:
        range_tot = np.append(range_tot, rslice)
        azi_tot = np.append(azi_tot, azislice)

    # Compute frequency map
    clutter_r, clutter_azi, freq = cmask_code.compute_frequency_map(rrange, azimuth, range_tot, azi_tot, nbfile)
    print("Clutter frequency map created.")

    # Some metadatas for the mask saved files.
    metakeys = ["site_name", "instrument_name", "author", "institution", "instrument_type"]
    metadata_out = dict()
    metadata_out['description'] = "Clutter mask"
    for mykey in metakeys:
        try:
            metadata_out[mykey] = radar.metadata[mykey]
        except KeyError:
            continue

    # Write mask to netcdf file.
    write_ncfile(outfilename_save, clutter_r, clutter_azi, metadata_out)
    print("Mask saved {}.".format(outfilename_save))

    if PLOT_FIG:
        try:
            plot_freq_map(outfilename_fig, rrange, freq)
            print("Figure saved {}.".format(outfilename_fig))
        except Exception:
            traceback.print_exc()
            pass

    return None


if __name__ == '__main__':
    """
    Global variables definition
    """
    DBZ_FIELD_NAME = "DBZ"
    RHOHV_FIELD_NAME = "RHOHV"
    INPUT_DIR = "/g/data2/rr5/vhl548/CPOL_level_1/2006/20060101/"
    OUTPUT_DIR = os.path.abspath("../saved_mask/")
    NCPU = 16
    PLOT_FIG = True

    warnings.simplefilter('ignore')
    main()
