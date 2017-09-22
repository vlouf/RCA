# Relative Calibration Adjusment technique

This is a ground radar calibration technique using the monitoring of ground clutters.

## Theory

The radar equation is:

$$ Z = 10 \log C + 20 \log r + 10 \log P_t$$

where $Z$ (dBZ) is the reflectivity factor, $r$ is the range, $P_t$ is the total power returned by the target, and $C$ is the so-called radar constant. To calibrate radars, we need to determine the value of $C$.

Ground echoes are generally caused by buildings, roads, or topographic structure near the radar. They are fixed targets $(\Delta r = 0)$, and their microphysics don't change over time $(\Delta P_t = 0)$, thus:

$$ \Delta Z_c = \Delta C $$

So any change in the ground clutter reflectivity is due to a change in the radar calibration.

## Usage

This is a 2-step technique: first create a clutter mask, then monitor ground clutter reflectivity.

### Step one: mask creation

Run (inside script directory) this line:

``` python RCA_step_one.py --input X --output X --rhohv X --dbz X --figure X --cpu X```

with:
>- ```--input```: Input directory for radar data (one day without precipitations is enough for creating the mask).
>- ```--output```: Ouput directory for saving the mask and/or the figure.
>- ```--rhohv```: Name of the cross-correlation ratio field (RHOHV by default).
>- ```--dbz```: Name of the uncorrected reflectivity (total power) field (DBZ by default).
>- ```--figure```: True/False if you want to plot the figure of the clutter frequency map (True by default).
>- ```--cpu```: Number of process used for reading input radar files (16 by default).

Only the `--input` argument is required.

This script will create a mask file (in output directory in netcdf format). This mask file is required for the second step.

### Step two: clutter monitoring

This step will use the mask from step 1 to extract the ground clutter reflectivity and compute the RCA value.
Run (inside script directory) this line:

``` python RCA_step_two.py --input X --output X --clutter X --dbz X --figure X --cpu X```

with:
>- ```--input```: Input directory for radar data (it can be an entire season).
>- ```--output```: Ouput directory for saving the mask and/or the figure.
>- ```--clutter```: Clutter mask file (the netcdf file created in step one).
>- ```--dbz```: Name of the uncorrected reflectivity (total power) field (DBZ by default).
>- ```--figure```: True/False if you want to plot the RCA time series (True by default).
>- ```--cpu```: Number of process used for reading input radar files (16 by default).

Only the `--input` and `--clutter` arguments are required.

## Requirements:

Required:
>- [Python 3][3] (tested with Python 3.5 and 3.6, probably - certainly - won't work with Python 2.)
>- Python ARM Radar Toolkit [(Py-ART) ][1]
>- Python Data Analysis Library [(pandas)][2]
>- netCDF4

Optionnal:
>- trmm_rsl library (for lassen files support)

The easiest way is to use [miniconda][3] to install the required python packages. Once miniconda is installed, type:

```conda install -c conda-forge arm_pyart trmm_rsl pandas numba```

This will install all required dependencies.

[1]: https://github.com/ARM-DOE/pyart
[2]: http://pandas.pydata.org/
[3]: https://conda.io/miniconda.html

## Note

There is a known bug in h5py (the HDF5 python library)  when reading an HDF5 ODIM file. This is corrected in the latest version of Py-ART (v1.9 still in development), but not yet in its released version (v1.8). An easy fix can be applied, just change the file odim_h5.py l.109 in /path/to/pyart/aux_io/odim_h5.py from: ```hfile = h5py.File(filename, 'r')``` to ```hfile = h5py.File(filename, 'r+', fclose_degree=h5py.h5f.CLOSE_DEFAULT) ```
