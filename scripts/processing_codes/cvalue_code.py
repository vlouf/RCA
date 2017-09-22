# Python standard library
import os
import datetime
import warnings
import traceback

from multiprocessing import Pool

# Other modules. Matplotlib must be imported first
import matplotlib
matplotlib.use('Agg')  # <- Reason why matplotlib is imported first.
import matplotlib.pyplot as pl
import pyart
import numpy as np
import pandas as pd


def compute_95th_percentile(the_data):
    """
    Compute the 95th percentile.

    Parameter:
    ==========
        the_data: np.array
            The clutter data.

    Returns:
    ========
        to_return: float
            95th percentile of the input data.
    """

    clut_series = pd.Series(the_data[~np.isnan(the_data)])
    to_return = clut_series.quantile(0.95)

    return to_return


def extract_clutter(r, azi, r_mask, th_mask, reflec, zdr):
    """
    Extract reflectivity data that has the same (r, azi) as the clutter mask.

    Parameters:
    ===========
        radar: struct
            Py-ART radar structure.
        r_mask: numpy.array(float)
            Clutter range.
        th_mask: numpy.array(float)
            Clutter azimuth.
        dbz_name: str
            Reflectivity field name.

    Returns:
    ========
        clut: numpy.array(float)
            Extracted clutter reflectivity.
    """
    reflec[reflec < 10] = np.NaN
    dr = r[1] - r[0]
    # Angle variation tolerance
    dazi = 0.5

    clut = np.array([])
    clut_zdr = np.array([])

    for the_r, the_azi in zip(r_mask, th_mask):
        pos_r = np.where((the_r >= r - dr) & (the_r < r + dr))[0]
        pos_azi = np.where((azi >= the_azi - dazi) &
                           (azi < the_azi + dazi))[0]

        if (len(pos_r) > 0) & (len(pos_azi) > 0):  # Non-empty array
            if len(pos_r) == len(pos_azi):
                rtmp = pos_r
                ttmp = pos_azi
            else:
                rtmp = pos_r[0]
                ttmp = pos_azi[0]
        else:
            continue

        try:
            value_rca = reflec[ttmp, rtmp]
            clut = np.append(clut, value_rca)
        except IndexError as err:
            print("Could not apply clutter mask.")
            traceback.print_exc()

        if zdr is not None:
            try:
                rca_zdr = zdr[ttmp, rtmp]
                clut_zdr = np.append(clut_zdr, rca_zdr)
            except IndexError as err:
                print("Could not apply clutter mask on ZDR.")
                traceback.print_exc()

    # Removing NAN values
    clut = clut[~np.isnan(clut)]
    clut_zdr = clut_zdr[~np.isnan(clut_zdr)]

    return clut, clut_zdr
