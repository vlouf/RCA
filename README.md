# Relative Calibration Adjusment technique

## Disclaimer

This dataset is supported by a funding from the U.S. Department of Energy as part of the Atmospheric Radiation Measurement (ARM) Climate Research Facility, an Office of Science user facility.

If you use this dataset to prepare a publication, please consider offering me (Valentin Louf) co-authorship and add the following line in the acknowledgments:

> This work has been supported by the U.S. Department of Energy Atmospheric Systems Research Program through the grant DE-SC0014063.

## Reference

- Louf, V., A. Protat, R. A. Warren, S. M. Collis, D. B. Wolff, S. Raunyiar, C. Jakob, and W. A. Petersen, 2018: An integrated approach to weather radar calibration and monitoring using ground clutter and satellite comparisons. J. Atmos. Ocean. Technol., JTECH-D-18-0007.1, doi:10.1175/JTECH-D-18-0007.1. [http://journals.ametsoc.org/doi/10.1175/JTECH-D-18-0007.1]

## Theory

This is a ground radar calibration technique using the monitoring of ground clutters.

The radar equation is:

![](http://latex.codecogs.com/gif.latex?Z%20%3D%2010%20%5Clog%20C%20&plus;%2020%20%5Clog%20r%20&plus;%2010%20%5Clog%20P_t)

where *Z* (dBZ) is the reflectivity factor, *r* is the range,  *Pt* is the total power returned by the target, and *C* is the so-called radar constant. To calibrate radars, we need to determine the value of *C*.

Ground echoes are generally caused by buildings, roads, or topographic structure near the radar. They are fixed targets ![](http://latex.codecogs.com/gif.latex?%28%5CDelta%20r%20%3D%200%29), and their microphysics don't change over time ![](http://latex.codecogs.com/gif.latex?%28%5CDelta%20P_t%20%3D%200%29), thus:

![](http://latex.codecogs.com/gif.latex?%5CDelta%20Z%20%3D%20%5CDelta%20C)
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
