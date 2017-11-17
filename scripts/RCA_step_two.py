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
    get_rain
    multproc_buffer_rca
    unpack_data
    main
"""
# Python standard library
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
import netCDF4
import numpy as np

# Custom modules.
from processing_codes import cvalue_code
from processing_codes import raijin_tools
from processing_codes.io import read_data, write_ncfile


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
    pl.plot(xdate, rca, '+')
    pl.plot([xdate.min(), xdate.max()], [rca.mean(), rca.mean()], 'r--')
    pl.xlim([xdate.min(), xdate.max()])
    pl.title("Monitoring of clutter reflectivity")
    pl.ylabel("CDF[clutter, 95%] (dBZ)")
    fig.autofmt_xdate()
    pl.savefig(outfilename_fig)
    pl.close()

    return None


def get_rain(r, azi, dbz, rhohv):
    """
    Compute a rough estimation of the rainfall rate at radar location (< 5 km)
    using an inverse Marshall-Palmer Z-R law.

    Parameters:
    ===========
    r: ndarray
        Range
    azi: ndarray
        Azimuth
    dbz: ndarray
        Reflectivity
    rhohv: ndarray
        Cross-correlation ratio

    Returns:
    ========
    rain_at_radar: float
        Estimated rainfall rate at radar site.
    """
    mydbz = dbz.copy()

    [R, T] = np.meshgrid(r, azi)
    # Keeping points at 5 km from radar & with RHOHV > 0.9
    mydbz[rhohv < 0.9] = np.NaN
    mydbz = mydbz[R < 5e3]
    rain = (10 ** (mydbz / 10) / 300) ** (2 / 3)
    rain_at_radar = np.nanmean(rain)

    return rain_at_radar


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
    volume_date, r, azi, reflec, zdr, rhohv = read_data(infile, DBZ_FIELD_NAME, ZDR_FIELD_NAME, RHOHV_FIELD_NAME)
    if rhohv is not None:
        rain_at_radar = get_rain(r, azi, reflec, rhohv)
    else:
        rain_at_radar = np.NaN

    try:
        ext_clut, clut_zdr = cvalue_code.extract_clutter(r, azi, range_permanent_echoes, azi_permanent_echoes, reflec, zdr)
    except Exception:
        print("Problem with this file:", os.path.basename(infile))
        traceback.print_exc()
        return volume_date, np.NaN, np.NaN, np.NaN

    rca = cvalue_code.compute_95th_percentile(ext_clut)
    if ZDR_FIELD_NAME is None:
        rca_zdr = np.NaN
    else:
        rca_zdr = cvalue_code.compute_95th_percentile(clut_zdr)

    return volume_date, rca, rca_zdr, rain_at_radar


def unpack_data(rslt):
    """
    Unpack data from multiprocessing and convert them to numpy array.

    Parameter:
    ==========
    rslt: list
        list of tuples returned by multiprocessing.

    Returns:
    ========
    xdate: ndarray
        Datetime array.
    rca: ndarray
        RCA values.
    rca_zdr: ndarray
        RCA for ZDR values.
    rain: ndarray
        Rainrate estimation at radar site.
    """
    # Unpack rslt
    xdate, rca, rca_zdr, rain = zip(*rslt)

    # Converting to numpy array
    xdate = np.array(xdate, dtype='datetime64[s]')
    rca = np.array(rca)
    rain = np.array(rain)
    # Sorting xdate (and rca) by chronological order.
    pos = np.argsort(xdate)
    xdate = xdate[pos]
    rca = rca[pos]
    rain = rain[pos]
    if ZDR_FIELD_NAME is not None:
        rca_zdr = np.array(rca_zdr)
        rca_zdr = rca_zdr[pos]

    return xdate, rca, rca_zdr, rain


def main():
    """
    1/ Generate file list.
    2/ Generate argument list for multiprocessing.
    3/ Invoke multiprocessing.
    4/ Unpack results from multiprocessing return and convert it to proper type.
    5/ Create output file name.
    6/ Save data.
    7/ (optionnal) plot figure.

    Parameters:
    ===========
    INPUT_DIR: str
        Radar data input directory.
    OUTPUT_DIR: str
        Clutter mask file.
    CLUTTER_RANGE: ndarray
        Clutter mask range.
    CLUTTER_AZIMUTH: ndarray
        Clutter mask azimuth.
    PLOT_FIG: bool
        Plot figure (True of False).
    NCPU: int
        Number of process
    """
    # Generating file list.
    flist = raijin_tools.get_files(INPUT_DIR)
    if len(flist) == 0:
        print("No file found.")
        return None
    print("Found %i files." % (len(flist)))

    # Create argument list for multiprocessing
    args_list = [(onefile, CLUTTER_RANGE, CLUTTER_AZIMUTH) for onefile in flist]

    # Calling the multiprocessing.
    with Pool(NCPU) as pool:
        rslt = pool.starmap(multproc_buffer_rca, args_list)
    print("Processing done for %i files. Unpacking data." % (len(rslt)))

    xdate, rca, rca_zdr, rain = unpack_data(rslt)

    # Output suffix str:
    st = xdate.min().tolist().isoformat()  # convert numpy datetime64 to datetime.datetime.
    ed = xdate.max().tolist().isoformat()
    # Data and figure output file names.
    outfilename = "RCA_{}_{}_to_{}.nc".format(INST_NAME, st, ed)
    outfilename = os.path.join(OUTPUT_DIR, outfilename)
    outfilename_fig = outfilename.replace(".nc", ".png")

    # Check if output file exists.
    if os.path.isfile(outfilename):
        print("Output file already exists {}. Doing nothing.".format(outfilename))
        return None

    gnrl_meta = dict()
    gnrl_meta['instrument_name'] = INST_NAME
    gnrl_meta['start_date'] = st
    gnrl_meta['end_date'] = ed
    write_ncfile(outfilename, xdate, rca, rca_zdr, rain, gnrl_meta)

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

    Global variables:
    =================
    INPUT_DIR: str
        Radar data input directory.
    OUTPUT_DIR: str
        Clutter mask file.
    CLUTTER_MASK_FILE: str
        Output directory.
    DBZ_FIELD_NAME: str
        Raw reflectivity (ZH) field name.
    ZDR_FIELD_NAME: str
        Differential reflectivity (ZDR) field name.
    PLOT_FIG: bool
        Plot figure (True of False).
    NCPU: int
        Number of process
    """
    # Argument parser.
    parser = argparse.ArgumentParser(description="RCA step 2: computing RCA results.")
    parser.add_argument('-i', '--input', dest='indir', default=None, type=str, help='Radar data input directory.')
    parser.add_argument('-c', '--clutter', dest='clutfile', default=None, type=str, help='Clutter mask file.')
    parser.add_argument('-o', '--output', dest='output', default=os.path.abspath("../saved_rca/"), type=str, help='Output directory.')
    parser.add_argument('-d', '--dbz', dest='dbz_name', default="DBZ", type=str, help='Raw reflectivity (ZH) field name.')
    parser.add_argument('-z', '--zdr', dest='zdr_name', default=None, type=str, help='Differential reflectivity (ZDR) field name.')
    parser.add_argument('-r', '--rhohv', dest='rhohv_name', default="RHOHV", type=str, help='Cross correlation ratio name.')
    parser.add_argument('-f', '--figure', dest='l_fig', default=False, type=bool, help='Plot figure (True of False).')
    parser.add_argument('-j', '--cpu', dest='ncpu', default=16, type=int, help='Number of process')

    # Global variables initialization.
    args = parser.parse_args()
    INPUT_DIR = args.indir
    OUTPUT_DIR = args.output
    CLUTTER_MASK_FILE = args.clutfile
    DBZ_FIELD_NAME = args.dbz_name
    ZDR_FIELD_NAME = args.zdr_name
    RHOHV_FIELD_NAME = args.rhohv_name
    PLOT_FIG = args.l_fig
    NCPU = args.ncpu

    # Checking global variables.
    if INPUT_DIR is None:
        parser.error("Need to provide an input directory.")

    # Check if input directory exists.
    if not os.path.isdir(INPUT_DIR):
        parser.error("Invalid input directory.")

    # Checking if the clutter mask is provided.
    if CLUTTER_MASK_FILE is None:
        parser.error("A clutter mask file is required.")

    # Checking if the output directory exists. Creating it otherwise.
    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
        print("Creating output directory.", OUTPUT_DIR)

    # Reading clutter mask.
    with netCDF4.Dataset(CLUTTER_MASK_FILE, "r") as ncid:
        try:
            CLUTTER_RANGE = ncid['range'][:]
            CLUTTER_AZIMUTH = ncid['azimuth'][:]
            INST_NAME = ncid.instrument_name
        except Exception:
            traceback.print_exc()
            print("Problem with the Clutter mask. Did you use RCA_step_one.py to create it?")
            sys.exit()

    if CLUTTER_RANGE is None:
        print("Clutter mask invalid.")
        sys.exit()

    warnings.simplefilter('ignore')
    main()
