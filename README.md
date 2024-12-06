# Synthetics

The original aim of this project is to provide an end-to-end pipeline to generate a set of synthetic LOFAR observations of a galaxy cluster, 
starting from a simulation cube of the cluster thermal emission. The output .MS file(s) contains the mock visibilities and,
unless specified otherwise, a column in which the visibilities have been corrupted by various effects, such as TEC, clock drifts, beam etc.
Images of the clean and corrupted visibilities are also automatically produced. 

## How to install

The code currently works through a singularity file which contains all the required packages. Only yt and LoSiTo (https://github.com/darafferty/losito)
needs to be installed separately, as they are not included in the container. If you wish to get the singularity image, please get in touch
with the developers.

## How to use

The preferred way is to use the singularity as described above. Enter in the singularity and move to your working directory.
Before starting, it is necessary to add the library directory to the PYTHONPATH:

`export PYTHONPATH=/path/to/Synthetics/synth_libs:$PYTHONPATH`

For LoSiTo to work correctly, one should also add:

`export PATH=/path/to/losito/bin:$PATH`

`export PYTHONPATH=/path/to/losito:$PYTHONPATH`

The code works through a simple command line which provides all the required information to produce the synthetic observations.
An example command could be:

`python /path/to/Synthetics/synthetics.py -p /path/to/simulation/cube -radec RA DEC`

These two inputs represent:

`-p`: this is the path to the simulation cube you want to start from to generate your mock observation.

`-radec`: these are the RA and DEC (in degrees) to use as phase center of your final .MS file.


These are the only required inputs, every other option you may add to the command line already has a default value and may be changed
if needed. These include:

`--name`: name of the output .MS file(s). Default is `simulated`. 

`--begin`: MJD time to start the mock observation from. 

`--duration`: Total observation time in hours. Default is 1 hour.

`--station`: Whether to use LBA or HBA. Default is LBA.

`--version`: Whether to use LOFAR1.0 or LOFAR2.0. Only LOFAR1.0 is currently tested.

`--freqrange`: The desired frequency range of the mock observation. Default is 42-66 MHz with LBA.

`--chanpersb`: Number of channels per subband. Default is 1.

`--chout`: Number of output channels when imaging. Default is 6.

`--imsize`: Image size. Default is 1024x1024.

`--nocorrupt`: If this option is given, no corruption is applied to the visibilities.

`--corruption_type`: Type of corruption to apply to visibilities. Options are: `tec, fr, clock, polmisalign, beam, 
noise, bandpass, all`. In the current version of the code, the user can specify the corruptions they wish to apply by simply separating them with spaces
(e.g. `--corruption_type tec fr`). If no `--corruption_type` is given as input, all types of corruptions are applied (unless `--nocorrupt` is specified).

`--recorrupt`: Use this option if you have already finished a run and you just want to change the type of corruption to apply to the 
visibilities.

`--skymodel_bdsf`: Use this option if you wish to use PyBDSF to get the sky model, instead of predicting from a .fits image. Necessary if one wish to apply
direction-dependent corruptions.

## Contributions

Most of the libraries used in this code are slightly modified versions of the libraries from `LiLF` (https://github.com/revoltek/LiLF). 
It also makes large use of `LoSiTo` (https://github.com/darafferty/losito) and `PyBDSF` (https://github.com/lofar-astron/PyBDSF). All credits go to their developers.


 