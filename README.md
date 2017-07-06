# Relative Calibration Adjusment technique

This is a ground radar calibration technique using the monitoring of ground clutters. By monitoring the reflectivity of ground clutters close to the radar, we are able to observe change in the radar calibration.

## Theory

The radar equation is:

$$ Z = 10 \log C + 20 \log r + 10 \log P_t$$

where $Z$ (dBZ) is the reflectivity factor, $r$ is the range, $P_t$ is the total power returned by the target, and $C$ is the so-called radar constant. To calibrate radars, we need to determine the value of $C$.

Ground echoes are generally caused by buildings, roads, or topographic structure near the radar. They are fixed targets $(\Delta r = 0)$, and their microphysics don't change over time $(\Delta P_t = 0)$, thus:

$$ \Delta Z_c = \Delta C $$

## Usage

This is a 2-step technique: first create a clutter mask, then monitor ground clutter reflectivity.

### Step one: mask creation

Run (inside script directory) this line:

``` python RCA_step_one.py --input --output --rhohv --dbz --figure --cpu```

with:
>- ```--input```: Input directory for radar data (one day without precipitations is enough for creating the mask).
>- ```--output```: Ouput directory for saving the mask and/or the figure.
>- ```--rhohv```: Name of the cross-correlation ratio field (RHOHV by default).
>- ```--dbz```: Name of the uncorrected reflectivity (total power) field (DBZ by default).
>- ```--figure```: True/False if you want to plot the figure of the clutter frequency map (True by default).
>- ```--cpu```: Number of process used for reading input radar files (16 by default).

This script will create a mask file (in output directory in netcdf format). This mask file is required for the second step.

### Step two: clutter monitoring

This step will use the mask from step 1 to extract the ground clutter reflectivity and compute the RCA value.
Run (inside script directory) this line:

``` python RCA_step_two.py --input --output --rhohv --dbz --figure --cpu```

with:
>- ```--input```: Input directory for radar data (it can be an entire season).
>- ```--output```: Ouput directory for saving the mask and/or the figure.
>- ```--clutter```: Clutter mask file (the netcdf file created in step one).
>- ```--dbz```: Name of the uncorrected reflectivity (total power) field (DBZ by default).
>- ```--figure```: True/False if you want to plot the RCA time series (True by default).
>- ```--cpu```: Number of process used for reading input radar files (16 by default).

## Requirements:

Mandatory:
>- Python 3 (tested with Python 3.5 and 3.6)
>- Py-ART
>- Pandas

Optionnal:
>- trmm_rsl library (for lassen files support)
