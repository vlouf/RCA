"""
Compute the RCA.

title: rca.py
author: Valentin Louf
email: valentin.louf@bom.gov.au
institution: Monash University and Bureau of Meteorology
date: 24/03/2020
"""
import gc

import pyart
import netCDF4
import numpy as np
import xarray as xr


def _read_radar(infile, refl_name):
    """
    Read input radar file

    Parameters:
    ===========
    radar_file_list: str
        List of radar files.
    refl_name: str
        Uncorrected reflectivity field name.

    Returns:
    ========
    radar: PyART.Radar
        Radar data.
    """
    try:
        radar = pyart.aux_io.read_odim_h5(infile, include_fields=[refl_name])
    except KeyError:
        gc.collect()  # Close file if stayed opened.
        radar = pyart.io.read(infile, include_fields=[refl_name])

    return radar


def extract_clutter(infile, mask_file, refl_name="total_power"):
    """
    Extract the clutter and compute the RCA value.
    
    Parameters:
    ===========
    infile: str
        Input radar file.
    mask_file: str
        Clutter mask dataset.
    refl_name: str
        Uncorrected reflectivity field name.
    
    Returns:
    ========
    dtime: np.datetime64
        Datetime of infile
    rca: float
        95th percentile of the clutter reflectivity.
    """
    # Radar data.
    radar = _read_radar(infile, refl_name)

    dtime = netCDF4.num2date(radar.time["data"][0], radar.time["units"])

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
    maskset = xr.open_dataset(mask_file)
    RC, AC = np.meshgrid(maskset.range, maskset.azimuth)

    npos = np.where(maskset.clutter_mask.values)
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
