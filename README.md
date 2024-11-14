# Synthetics

This project aims to provide an end-to-end pipeline to generate a set of synthetic LOFAR observations of a galaxy cluster, 
starting from a simulation cube of the cluster thermal emission. The output .MS file(s) contains the mock visibilities and,
if needed, a column in which the visibilities have been corrupted by various effects, such as TEC, clock drifts, beam etc.
Images of the clean and corrupted visibilities are also automatically produced.

## How to install

The code currently works through a singularity file which contains all the required packages. Only yt and LoSiTo (https://github.com/darafferty/losito)
needs to be installed separately, as they are not included in the container.

## How to use

The preferred way is to use the singularity as described above, enter in the singularity and move to your working directory.
Before starting, it is necessary to add the library directory to the PYTHONPATH:

`export PYTHONPATH=/path/to/Synthetics/synth_libs:$PYTHONPATH`

You may also want to add the Synthetics directory to PATH.
The code works through a simple command line which provides all the required information to produce the synthetic observations.
An example command could be:

`python /path/to/Synthetics/synthetics.py -synth /path/to/Synthetics -p /path/to/simulation/cube -radec RA DEC`

These three inputs represent:

`-synth`: this is simply the path to the Synthetics directory you cloned from GitHub.

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

`--corruption_type `: Type of corruption to apply to visibilities. Options are: `tec, tec_fr, tec_fr_clock, polmisalign, beamcorrupt, 
noise, all`. In the current version of the code, these corruptions are in 'ascending' order: each corruption also includes all the 
previous ones (e.g. polmisalign will include TEC, Faraday roation and clock drifts). Default is `all`.

`--recorrupt`: Use this option if you have already finished a run and you just want to change the type of corruption to apply to the 
visibilities.


## Contributions

Most of the libraries used in this code are slightly modified versions of the libraries from `LiLF` (https://github.com/revoltek/LiLF). 
It also makes great use of `LoSiTo` (https://github.com/darafferty/losito). All credits go to their developers.


 