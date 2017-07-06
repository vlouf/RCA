import pyart
import numpy as np

from numba import jit


@jit(nopython=True)
def _jit_find_clut_pos(rrange, theta, range_list, azi_list):
    # Initializing frequency array
    frequency = np.zeros((len(theta), len(rrange)))
    # Running through lists of non-meteorological echoes
    for myr, mya in zip(range_list, azi_list):
        # Finding their positions
        apos = np.argmin(np.abs(theta - mya))
        rpos = np.argmin(np.abs(rrange - myr))

        frequency[apos, rpos] += 1
    return frequency


def get_clutter_position(radar,
                         dbz_name="DBZ",
                         rhohv_name="RHOHV",
                         refl_thrld=45,
                         rhohv_thrld=0.6,
                         maxrange=10e3):
    """
    Get (range, azimuth) position of non meteorological echoes.

    Parameters:
    ===========
        radar: struct
            Py-ART radar structure.
        dbz_name: str
            Reflectivity field name.
        rhohv_name: str
            Cross-correlation field name.
        refl_thrld: float
            Minimum reflectivity threshold (in dBZ)
        rhohv_thrld: float
            Maximum cross correlation threshold below which echoes
            are not from a meteorological source.
        maxrange: int
            Maximum range (in meters)

    Returns:
    ========
        r_clutter: np.array
            Range position of non-meteorological echoes.
        azi_clutter: np.array
            Azimuth position of non-meteorological echoes.
    """
    # Extract first elevation only
    rslice = radar.get_slice(0)
    # Extract range/azimuth
    r = radar.range['data'].astype(int)
    azi = radar.azimuth['data'][rslice]
    # Check azimuth field.
    if len(azi) <= 60:
        print("Invalid azimuth field")
        return None

    # Get reflectivity and RHOHV
    total_power = radar.fields[dbz_name]['data'][rslice].filled(np.NaN)
    cross_correlation_ratio = radar.fields[rhohv_name]['data'][rslice].filled(np.NaN)

    # Removing every echoes that are above RHOHV threshold, below DBZ threshold and above maximum range.
    clut = total_power
    clut[(cross_correlation_ratio > rhohv_thrld) | (total_power < refl_thrld)] = np.NaN
    clut[:, r > maxrange] = np.NaN

    # Position of remaining echoes
    posa, posr = np.where(~np.isnan(clut))

    r_clutter = r[posr]
    azi_clutter = azi[posa]

    return r_clutter, azi_clutter


def compute_frequency_map(rrange,
                          azimuth,
                          range_list,
                          azi_list,
                          nb_files,
                          freq_thrld=80):
    """
    Compute the frequency map of non-meteorological echoes

    Parameters:
    ===========
        rrange: np.array
            Radar range
        azimuth: np.array
            Radar azimuth
        range_list: np.array
            Range position of non-meteorological echoes
        azi_list: np.array
            Azimuth position of non-meteorological echoes
        nb_files: int
            Total number of files for which the range, azimuth of
            non-meteorological echoes as been retrieved.
        freq_thrld: int
            Frequency threshold in % above which echoes are
            considered permanent.

    Returns:
    ========
        clutter_range: np.array
            Estimated permanent clutter range
        clutter_azimuth: np.array
            Estimated permanent clutter azimuth
        frequency: np.array
            Frequency map of clutter echoes.
    """
    # Azimuthal width
    # Because radar azimuth don't always (never) starts at 0, we create
    # a theta angle array that starts at 0 with the same step as azimuth.
    da = azimuth[1] - azimuth[0]
    da = np.round(da, 1)  # Round to 1 figure after the decimal.

    # Generating position arrays
    theta = np.arange(0, 360, da)
    [TH, R] = np.meshgrid(theta, rrange, indexing='ij')

    frequency = _jit_find_clut_pos(rrange, theta, range_list, azi_list)

    # It can happen when da is not an integer.
    frequency[frequency > nb_files] = nb_files
    # In percent
    frequency = 100*frequency/nb_files

    if np.max(frequency) < 90:
        print("Not enough clutter. Try another day")
        return None

    # Finding permanent echoes azimuth and range position
    azpos, rapos = np.where(frequency > freq_thrld)
    clutter_range = rrange[rapos]
    clutter_azimuth = theta[azpos]

    return clutter_range, clutter_azimuth, frequency
