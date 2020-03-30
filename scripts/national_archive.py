'''
Radar calibration monitoring using ground clutter. Processing the Australian
National archive.

@creator: Valentin Louf <valentin.louf@bom.gov.au>
@institution: Monash University and Bureau of Meteorology
@date: 19/03/2020

    buffer
    check_rid
    extract_zip
    get_radar_archive_file
    mkdir
    remove
    savedata
    main
'''
import gc
import os
import sys
import glob
import time
import zipfile
import argparse
import datetime
import warnings
import traceback

import crayons
import numpy as np
import pandas as pd
import dask.bag as db

import cluttercal


def buffer(infile, cmask):
    '''
    Buffer function to catch and kill errors about missing Sun hit.

    Parameters:
    ===========
    infile: str
        Input radar file.

    Returns:
    ========
    rslt: pd.DataFrame
        Pandas dataframe with the results from the solar calibration code.
    '''
    try:
        dtime, rca = cluttercal.extract_clutter(infile, cmask, refl_name='total_power')
    except ValueError:
        return None
    except Exception:
        print(infile)
        traceback.print_exc()
        return None
    
    return dtime, rca


def check_rid():
    '''
    Check if the Radar ID provided exists.
    '''
    indir = f'/g/data/rq0/odim_archive/odim_pvol/{RID}'
    return os.path.exists(indir)


def extract_zip(inzip, path):
    '''
    Extract content of a zipfile inside a given directory.

    Parameters:
    ===========
    inzip: str
        Input zip file.
    path: str
        Output path.

    Returns:
    ========
    namelist: List
        List of files extracted from  the zip.
    '''
    with zipfile.ZipFile(inzip) as zid:
        zid.extractall(path=path)
        namelist = [os.path.join(path, f) for f in zid.namelist()]
    return namelist


def get_radar_archive_file(date):
    '''
    Return the archive containing the radar file for a given radar ID and a
    given date.

    Parameters:
    ===========
    date: datetime
        Date.

    Returns:
    ========
    file: str
        Radar archive if it exists at the given date.
    '''
    datestr = date.strftime('%Y%m%d')
    file = f"/g/data/rq0/odim_archive/odim_pvol/{RID}/{date.year}/vol/{RID}_{datestr}.pvol.zip"
    if not os.path.exists(file):
        return None

    return file


def mkdir(path):
    '''
    Create the DIRECTORY(ies), if they do not already exist
    '''
    try:
        os.mkdir(path)
    except FileExistsError:
        pass

    return None


def remove(flist):
    '''
    Remove file if it exists.
    '''
    flist = [f for f in flist if f is not None]
    for f in flist:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass
    return None


def savedata(df, date, path):
    '''
    Save the output data into a CSV file compatible with pandas.

    Parameters:
    ===========
    df: pd.Dataframe
        RCA timeserie to be saved.
    date:
        Date of processing.
    path: str
        Output directory.
    '''
    datestr = date.strftime('%Y%m%d')

    path = os.path.join(path, 'rca')
    mkdir(path)
    path = os.path.join(path, RID)
    mkdir(path)
    path = os.path.join(path, str(date.year))
    mkdir(path)

    outfilename = os.path.join(path, f'rca.{RID}.{datestr}.csv')
    df.to_csv(outfilename)
    print(crayons.green(f'{len(df)} solar hits on {datestr}.'))
    print(crayons.green(f'Results saved in {outfilename}.'))

    return None


def gen_cmask(radar_file_list, date):
    '''
    Generate the clutter mask for a given day and save the clutter mask as a
    netCDF.

    Parameters:
    ===========
    radar_file_list: list
        List radar files for the given date.
    date: datetime
        Date.

    Returns:
    ========
    file_prefix: str
        Prefix given to filename.
    outpath: str
        Output directory for the clutter masks.
    '''
    file_prefix = f'{RID}_'
    datestr = date.strftime('%Y%m%d')

    outpath = os.path.join(OUTPATH, 'cmasks', f'{RID}')
    mkdir(outpath)
    outputfile = os.path.join(outpath, file_prefix + f'{datestr}.nc')

    if os.path.isfile(outputfile):
        print('Clutter masks already exists. Doing nothing.')
        return file_prefix

    cmask = cluttercal.clutter_mask(radar_file_list,
                                    refl_name="total_power",
                                    refl_threshold=50,
                                    max_range=20e3,
                                    freq_threshold=50,
                                    use_dask=True)
    cmask.to_netcdf(outputfile)

    return file_prefix, outpath


def main(date_range):
    for date in date_range:
        zipfile = get_radar_archive_file(date)
        if zipfile is None:
            print(crayons.red(f'No file found for date {date}.'))
            continue

        namelist = extract_zip(zipfile, path=ZIPDIR)
        print(crayons.yellow(f'{len(namelist)} files to process for {date}.'))

        # Generate clutter masks.
        prefix, outpath = gen_cmask(namelist, date)

        try:
            cmask = cluttercal.composite_mask(date, timedelta=7, indir=outpath, prefix=prefix)
        except ValueError:
            cmask = cluttercal.single_mask(date, indir=outpath, prefix=prefix)

        arglist = [(f, cmask) for f in namelist]

        bag = db.from_sequence(arglist).starmap(buffer)
        rslt = bag.compute()
        if rslt is None:
            print(crayons.red(f'No results for date {date}.'))
            remove(namelist)
            continue

        rslt = [r for r in rslt if r is not None]
        if len(rslt) != 0:
            ttmp, rtmp = zip(*rslt)
            rca = np.array(rtmp)
            dtime = np.array(ttmp, dtype='datetime64')

            if len(rca) == 0:
                print(crayons.red(f'No results for radar {RID}.'))
            else:
                df = pd.DataFrame({'rca': rca}, index=dtime)
                savedata(df, date, path=OUTPATH)
        else:
            print(crayons.red(f'No results for date {date}.'))

        # Removing unzipped files, collecting memory garbage.
        remove(namelist)
        del bag
        gc.collect()

    return None


if __name__ == "__main__":
    parser_description = "Solar calibration of radar in the National radar archive."
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        "-r",
        "--rid",
        dest="rid",
        type=int,
        required=True,
        help="Radar ID number.")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default='/scratch/kl02/vhl548/rca_output/',
        type=str,
        help="Output directory")
    parser.add_argument(
        '-s',
        '--start-date',
        dest='start_date',
        default=None,
        type=str,
        help='Starting date.',
        required=True)
    parser.add_argument(
        '-e',
        '--end-date',
        dest='end_date',
        default=None,
        type=str,
        help='Ending date.',
        required=True)

    args = parser.parse_args()
    RID = f"{args.rid:02}"
    START_DATE = args.start_date
    END_DATE = args.end_date
    OUTPATH = args.output
    ZIPDIR = '/scratch/kl02/vhl548/unzipdir/'

    if not check_rid():
        parser.error('Invalid Radar ID.')
        sys.exit()

    try:
        start = datetime.datetime.strptime(START_DATE, "%Y%m%d")
        end = datetime.datetime.strptime(END_DATE, "%Y%m%d")
        if start > end:
            parser.error('End date older than start date.')
        date_range = [start + datetime.timedelta(days=x) for x in range(0, (end - start).days + 1, )]
    except ValueError:
        parser.error('Invalid dates.')
        sys.exit()

    print(crayons.green(f'Processing sun calibration for radar {RID}.'))
    print(crayons.green(f'Between {START_DATE} and {END_DATE}.'))
    print(crayons.green(f'Data will be saved in {OUTPATH}.'))

    tick = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        main(date_range)
    tock = time.time()
    print(crayons.magenta(f'Process finished in {tock - tick:0.2}s.'))
