# Relative Calibration Adjusment technique

This is a ground radar calibration technique using the monitoring of ground clutters.

## Theory

The radar equation is:

<img src="http://www.sciweavers.org/tex2img.php?eq=%20Z%20%3D%2010%20%5Clog%20C%20%2B%2020%20%5Clog%20r%20%2B%2010%20%5Clog%20P_t&bc=White&fc=Black&im=jpg&fs=12&ff=fourier&edit=0" align="center" border="0" alt=" Z = 10 \log C + 20 \log r + 10 \log P_t" width="219" height="18" />

where *Z* (dBZ) is the reflectivity factor, *r* is the range,  *Pt* is the total power returned by the target, and *C* is the so-called radar constant. To calibrate radars, we need to determine the value of *C*.

Ground echoes are generally caused by buildings, roads, or topographic structure near the radar. They are fixed targets <img src="http://www.sciweavers.org/tex2img.php?eq=%28%5CDelta%20r%20%3D%200%29&bc=White&fc=Black&im=jpg&fs=12&ff=fourier&edit=0" align="center" border="0" alt="(\Delta r = 0)" width="58" height="17" />, and their microphysics don't change over time <img src="http://www.sciweavers.org/tex2img.php?eq=%28%5CDelta%20P_t%20%3D%200%29&bc=White&fc=Black&im=jpg&fs=12&ff=fourier&edit=0" align="center" border="0" alt="(\Delta P_t = 0)" width="65" height="17" />, thus:

<img src="http://www.sciweavers.org/tex2img.php?eq=%20%5CDelta%20Z_c%20%3D%20%5CDelta%20C&bc=White&fc=Black&im=jpg&fs=12&ff=fourier&edit=0" align="center" border="0" alt=" \Delta Z_c = \Delta C" width="71" height="17" />

So any change in the ground clutter reflectivity is due to a change in the radar calibration.

## Usage

This is a 2-step technique: first create a clutter mask, then monitor ground clutter reflectivity.

### Step one: mask creation

Run (inside script directory) this line:

``` python RCA_step_one.py --input X --output X --rhohv X --dbz X --figure X --cpu X```

with:
>- ```--input```: Input directory for radar data (one day without precipitations is enough for creating the mask).
>- ```--output```: Ouput directory for saving the mask and/or the figure.
>- ```--rhohv```: Name of the cross-correlation ratio field (RHOHV by default). The cross-correlation ratio is optional.
>- ```--dbz```: Name of the uncorrected reflectivity (total power) field (DBZ by default).
>- ```--figure```: True/False if you want to plot the figure of the clutter frequency map (True by default).
>- ```--cpu```: Number of process used for reading input radar files (16 by default).

Only the `--input` argument is required.

This script will create a mask file (in output directory in [netcdf][4] format). This mask file is required for the second step.

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
[4]: https://www.unidata.ucar.edu/software/netcdf/
