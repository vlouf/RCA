"""
Second part of the ground-clutter-monitoring calibration technique. This part
monitors the ground clutter reflectivity.

@title: RCA_step_two
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Bureau of Meteorology/Monash University
@date: 5/07/2017
@version: 0.99

.. autosummary::
    :toctree: generated/

    plot_figure
    write_ncfile
    multproc_buffer_rca
    main
"""
# Python standard library
import gc
import os
import sys
import argparse
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
from processing_codes import cvalue_code
from processing_codes.raijin_tools import get_files


def plot_figure(outfilename_fig, xdate, rca):
    """
    Plot the clutter frequency map.

    Parameters:
    ===========
        outfilename_fig: str
            Output file name for figure.
        xdate: np.array
            Datetime
        rca: np.array
            RCA value.
    """
    fig = pl.figure()
    pl.plot(xdate, rca)
    pl.plot([xdate.min(), xdate.max()], [rca.mean(), rca.mean()], 'r--')
    pl.xlim([xdate.min(), xdate.max()])
    pl.title("Monitoring of clutter reflectivity")
    pl.ylabel("CDF[clutter, 95%] (dBZ)")
    fig.autofmt_xdate()
    pl.savefig(outfilename_fig, dpi=150)
    pl.close()

    return None


def write_ncfile(outfilename, xdate, rca, gnrl_meta):
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

        # Assign values.
        ncr[:] = time
        nca[:] = rca

        # Set units.
        ncr.units = time_units
        nca.units = "degrees"

        # Set main metadata
        for mykey in gnrl_meta.keys():
            rootgrp.setncattr_string(mykey, gnrl_meta[mykey])

    return None


def multproc_buffer_rca(infile, range_permanent_echoes, azi_permanent_echoes):
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
        volume_date: datetime
            Datetime for input volume
        rca: float
            CDF[95%] of ground clutter reflectivity
    """
    try:
        radar = pyart.io.read(infile)
    except Exception:
        print("Could not read input file", os.path.basename(infile))
        return None

    volume_date = netCDF4.num2date(radar.time['data'][0], radar.time['units'])

    try:
        ext_clut = cvalue_code.extract_clutter(radar, range_permanent_echoes, azi_permanent_echoes, DBZ_FIELD_NAME)
    except Exception:
        print("Problem with this file:", os.path.basename(infile))
        traceback.print_exc()
        return None

    rca = cvalue_code.compute_95th_percentile(ext_clut)
    gc.collect()  # Free unused memory

    return volume_date, rca


def main():
    flist = get_files(INPUT_DIR)

    # Create argument list for multiprocessing
    args_list = []
    for onefile in flist:
        tmp = (onefile, CLUTTER_RANGE, CLUTTER_AZIMUTH)
        args_list.append(tmp)

    with Pool(NCPU) as pool:
        rslt = pool.starmap(multproc_buffer_rca, args_list)

    # Unpack rslt
    rca = []
    xdate = []
    for tmptime, tmprca in rslt:
        xdate.append(tmptime)
        rca.append(tmprca)

    # Sorting xdate (and rca) by chronological order.
    rca = np.array(rca)
    xdate = np.array(xdate, dtype='datetime64[s]')
    pos = np.argsort(xdate)
    xdate = xdate[pos]
    rca = rca[pos]

    # Output suffix str:
    st = xdate.min().tolist().isoformat()
    ed = xdate.max().tolist().isoformat()
    # Data output file name.
    outfilename = "RCA_{}_{}_to_{}.nc".format(INST_NAME, st, ed)
    outfilename = os.path.join(OUTPUT_DIR, outfilename)
    # Figure output file name.
    outfilename_fig = "RCA_{}_{}_to_{}.png".format(INST_NAME, st, ed)
    outfilename_fig = os.path.join(OUTPUT_DIR, outfilename_fig)

    # Check if output file exists.
    if os.path.isfile(outfilename):
        print("Output file already exists {}.".format(outfilename))
        print("Doing nothing.")
        return None

    gnrl_meta = dict()
    gnrl_meta['instrument_name'] = INST_NAME
    gnrl_meta['start_date'] = st
    gnrl_meta['end_date'] = ed
    write_ncfile(outfilename, xdate, rca, gnrl_meta)

    if PLOT_FIG:
        try:
            plot_figure(outfilename_fig, xdate, rca)
            print("Figure saved {}.".format(outfilename_fig))
        except Exception:
            traceback.print_exc()
            pass

    return None


if __name__ == '__main__':
    """
    Global variables definition
    """
    # Argument parser.
    parser = argparse.ArgumentParser(description="RCA step 1: creation of the clutter mask.")
    parser.add_argument('-i', '--input', dest='indir', default=None, type=str, help='Radar data input directory.')
    parser.add_argument('-c', '--clutter', dest='clutfile', default=None, type=str, help='Clutter mask file.')
    parser.add_argument('-o', '--output', dest='output', default=os.path.abspath("../saved_rca/"), type=str, help='Output directory.')
    parser.add_argument('-d', '--dbz', dest='dbz_name', default="DBZ", type=str, help='Raw reflectivity (total power) field name.')
    parser.add_argument('-f', '--figure', dest='l_fig', default=True, type=bool, help='Plot figure (True of False).')
    parser.add_argument('-j', '--cpu', dest='ncpu', default=16, type=int, help='Number of process')

    # Global variables initialization.
    args = parser.parse_args()
    INPUT_DIR = args.indir
    OUTPUT_DIR = args.output
    CLUTTER_MASK_FILE = args.clutfile
    DBZ_FIELD_NAME = args.dbz_name
    PLOT_FIG = args.l_fig
    NCPU = args.ncpu

    # Checking global variables.
    if INPUT_DIR is None:
        parser.error("Need to provide an input directory.")

    if not os.path.isdir(INPUT_DIR):
        parser.error("Invalid input directory.")

    if CLUTTER_MASK_FILE is None:
        parser.error("A clutter mask file is required.")

    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
        print("Creating output directory.", OUTPUT_DIR)

    with netCDF4.Dataset(CLUTTER_MASK_FILE, "r") as ncid:
        CLUTTER_RANGE = ncid['range'][:]
        CLUTTER_AZIMUTH = ncid['azimuth'][:]
        INST_NAME = ncid.instrument_name

    if CLUTTER_RANGE is None:
        print("Clutter mask invalid.")
        sys.exit()

    warnings.simplefilter('ignore')
    main()
