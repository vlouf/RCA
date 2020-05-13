"""
Compute the RCA.

title: rca.py
author: Valentin Louf
email: valentin.louf@bom.gov.au
institution: Monash University and Bureau of Meteorology
date: 30/03/2020
"""
import gc
import os

import pyart
import cftime
import numpy as np
import pandas as pd
import xarray as xr

from .cluttermask import _read_radar


def composite_mask(date, timedelta=7, indir="compomask", prefix="cpol_cmask_", freq_thrld=0.9):
    """
    Generate composite clutter mask.

    Parameters:
    -----------
    date: Timestamp
        Date of processing
    timedelta: int
        Time delta for the composite.
    indir: str
        Where clutter mask are stored
    prefix: str
        What is the clutter mask file prefix.
    freq_thrld: float
        Frequency threshold (0 < f < 1).

    Returns:
    --------
    mask: ndarray
        Clutter mask
    """

    def get_mask_list(date, timedelta, indir, prefix):
        drange = pd.date_range(date - pd.Timedelta(f"{timedelta}D"), date)
        flist = []
        for day in drange:
            file = os.path.join(indir, prefix + "{}.nc".format(day.strftime("%Y%m%d")))
            if os.path.isfile(file):
                flist.append(file)
        return flist

    flist = get_mask_list(date, timedelta, indir, prefix)
    if len(flist) == 1:
        composite = xr.open_dataset(flist[0]).clutter_mask.values
    else:
        cmask = [xr.open_dataset(f).clutter_mask.values for f in flist]
        cmaskarr = np.concatenate(cmask, axis=np.newaxis).reshape((len(flist), 360, 20))
        compo_freq = cmaskarr.sum(axis=0) / len(flist)
        composite = compo_freq > freq_thrld
        if np.sum(composite) == 0:
            print(f'BAD COMPOSITE FOR {date}')
            composite = compo_freq != 0

    return composite


def single_mask(date, indir="compomask", prefix="cpol_cmask_"):
    """
    Generate clutter mask.

    Parameters:
    -----------
    date: Timestamp
        Date of processing
    indir: str
        Where clutter mask are stored
    prefix: str
        What is the clutter mask file prefix.

    Returns:
    --------
    mask: ndarray
        Clutter mask
    """
    file = os.path.join(indir, prefix + "{}.nc".format(date.strftime("%Y%m%d")))
    if not os.path.isfile(file):
        file = sorted(glob.glob(os.path.join(indir, '*.nc')))[-1]

    print(f'Using emergency mask {file}')
    cmask = xr.open_dataset(file).clutter_mask.values
    return cmask


def extract_clutter(infile, clutter_mask, refl_name="total_power"):
    """
    Extract the clutter and compute the RCA value.

    Parameters:
    -----------
    infile: str
        Input radar file.
    clutter_mask: numpy.array(float)
        Clutter mask (360 deg x 20 km)
    refl_name: str
        Uncorrected reflectivity field name.

    Returns:
    --------
    dtime: np.datetime64
        Datetime of infile
    rca: float
        95th percentile of the clutter reflectivity.
    """
    # Radar data.
    radar = _read_radar(infile, refl_name)

    dtime = cftime.num2pydate(radar.time["data"][0], radar.time["units"])

    sl = radar.get_slice(0)
    r = radar.range["data"]
    azi = radar.azimuth["data"][sl]
    try:
        refl = radar.fields[refl_name]["data"][sl][:, r < 20e3].filled(np.NaN)
    except AttributeError:
        refl = radar.fields[refl_name]["data"][sl][:, r < 20e3]
    zclutter = np.zeros_like(refl) + np.NaN

    r = r[r < 20e3]
    R, A = np.meshgrid(r, azi)
    R = (R // 1000).astype(int)
    A = (np.round(A) % 360).astype(int)

    # Mask.
    RC, AC = np.meshgrid(np.arange(20), np.arange(360))

    npos = np.where(clutter_mask)
    for ir, ia in zip(RC[npos], AC[npos]):
        pos = (R == ir) & (A == ia)
        zclutter[pos] = refl[pos]

    try:
        zclutter = zclutter[~np.isnan(zclutter)]
        rca = np.percentile(zclutter, 95)
    except IndexError:
        # Empty array full of NaN.
        raise ValueError("All clutter is NaN.")

    del radar
    return dtime, rca
