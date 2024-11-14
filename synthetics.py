import sys, os, glob, argparse
import numpy as np
from synth_libs import lib_ms, lib_img, lib_log, lib_util
import warnings
import astropy
from astropy.cosmology import FlatLambdaCDM
import yt
from yt.units import gauss, g, cm, Hz, W, m, Jy, erg, s, G
from astropy.io import fits
import math

logger_obj = lib_log.Logger('synthetics')
logger = lib_log.logger
sch = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('synthetics.walker')
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale = '3asec'

##################################################################

def _synchrotron_emissivity(field, data):
    #n_e = data["gas", "density"]/9.11e-28 # Converting density into numeric density #TODO this is gas density, we need the relativistic electron density
    n_e = 1e-5 * cm**-3
    p = 2.5  # Power-law index (2.5 is quite typical...)
    B_perp = abs(data["gas", "magnetic_field_y"]*1e-6) #We need this 1e-6 factor because for some reason it's missing when read from within the function...
    emissivity = n_e * (B_perp** ((p+1)/2)) / G**(7/4) * erg * s**-1 * cm**-2 * cm**3 #We add the frequency dependency later because it's not in the simulation
    return emissivity

def convert_coordinates(lat_deg, lon_deg):
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    return lat_rad, lon_rad

##################################################################

parser = argparse.ArgumentParser(description='Create synthetic .MS file which contains the visibilities of the synchrotron emission derived from a simulation cube, as seen by the LOFAR telescope.')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for the simulation cube.')
parser.add_argument('-synth', '--synth', dest='synthpath', action='store', default='', type=str, help='Path to the Synthetics directory cloned from GitHub. Example: /path/to/Synthetics/')
parser.add_argument('-radec', '--radec', dest='radec', nargs=2, type=float, default=None, help='RA/DEC of the phase centre of the generated .MS file.')
parser.add_argument('--name', dest='name', type=str, default='simulated', help='Name of the generated .MS file.')
parser.add_argument('--begin', dest='begin', type=float, default='60310', help='MJD starting time of the generated .MS file.')
parser.add_argument('--duration', dest='duration', type=float, default=1, help='Duration of the generated .MS file in hours. Default is 1 hour.')
parser.add_argument('--station', dest='station', type=str, default='LBA', help='Whether to use LBA or HBA. Default is LBA.')
parser.add_argument('--version', dest='version', type=int, default=1, help='Whether to use LOFAR or LOFAR 2.0. Default is LOFAR. Use 2 for LOFAR 2.0')
parser.add_argument('--freqrange', dest='freqrange', nargs=2, type=float, default=(42e6, 66e6), help='Frequency range of the generated .MS file in Hz.')
parser.add_argument('--chanpersb', dest='chanpersb', type=int, default=1, help='Number of channels per sub-band. Default is 1.')
parser.add_argument('--chout', dest='chout', type=int, default=6, help='Number of output channel images in WSClean. Default is 6.')
parser.add_argument('--imsize', dest='imsize', type=int, default=1024, help='Pixel size of images. Default is 1024.')
parser.add_argument('--nocorrupt', dest='nocorrupt', action='store_true', help='Whether to corrupt the dataset with ionospheric delays.')
parser.add_argument('--corruption_type', dest='corrtype', action='store', default='all', type=str, help='Type of corruption to apply to the dataset. Can be set to "tec", "tec_fr", "tec_fr_clock", "polmisalign", "beamcorrupt", "noise" or "all", in ascending order (i.e. the latter includes all the formers. Default is "all".')
parser.add_argument('--recorrupt', dest='recorrupt', action='store_true', help='Use this if you just want to change the type of corruption to apply, without re-running everything else.')

args = parser.parse_args()
pathdir = args.path
deg_coords = args.radec
synth = args.synthpath
name = args.name
start = args.begin
duration = args.duration
station = args.station
version = args.version
freqrange = args.freqrange
chanpersb = args.chanpersb
chout = args.chout
imsize = args.imsize
nocorrupt = args.nocorrupt
corrtype = args.corrtype
recorrupt = args.recorrupt

#Since we allow coordinates in degrees, we need to convert to radians for synthms
coords = convert_coordinates(deg_coords[0], deg_coords[1])

if corrtype not in ["tec", "tec_fr", "tec_fr_clock", "polmisalign", "beamcorrupt", "noise", "all"]:
    logger.error('Unknown corruption type. Possible values are tec, tec_fr, tec_fr_clock, polmisalign, beamcorrupt, noise or all.')
    sys.exit()

if not pathdir:
    logger.error('Provide a path (-p) where to look for the simulation cube.')
    sys.exit()

if not coords:
    logger.error('Provide RA and DEC in degrees of the phase centre of the .MS file you want to generate.')
    sys.exit()

if not synth:
    logger.error('Provide the path to the Synthetics directory cloned from GitHub.')
    sys.exit()

if corrtype == 'tec':
    parset = f'{synth}/parsets/tec.parset'
elif corrtype == 'tec_fr':
    parset = f'{synth}/parsets/tec_fr.parset'
elif corrtype == 'tec_fr_clock':
    parset = f'{synth}/parsets/tec_fr_clock.parset'
elif corrtype == 'polmisalign':
    parset = f'{synth}/parsets/polmisalign.parset'
elif corrtype == 'beamcorrupt':
    parset = f'{synth}/parsets/beamcorrupt.parset'
elif corrtype == 'noise':
    parset = f'{synth}/parsets/noise.parset'
elif corrtype == 'all':
    parset = f'{synth}/parsets/bandpass.parset'

with w.if_todo('cleaning'):
    logger.info('Preparing the environment...')
    lib_util.check_rm('data')
    lib_util.check_rm('images')
    os.makedirs('data')
    os.makedirs('images')
    lib_util.check_rm('models')
    os.makedirs('models')

with w.if_todo('generate_MS'):
    logger.info('Generating empty .MS files...')
    sch.add(f'synthms --name data/{name} --start {start} --tobs {duration} --ra {coords[0]} --dec {coords[1]} --station {station} '
          f'--lofarversion {version} --minfreq {freqrange[0]} --maxfreq {freqrange[1]} --chanpersb {chanpersb}', log='generateMS.log', commandType='general', processors='max')
    sch.run(check=True)

logger.info('Opening empty MS files...')
MSs_empty = lib_ms.AllMSs(glob.glob(f'data/{name}*.MS'), sch, check_flags=False)

with w.if_todo('empty_image'):
    logger.info('Producing empty image from initial dataset...')
    sch.add(f'wsclean -size {imsize} {imsize} -name images/empty -scale {scale} -data-column DATA -weight briggs -1 -circular-beam -niter 100000 -no-update-model-required -mgain 0.6'
          f' -baseline-averaging 10 -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
          log='wsclean-empty.log', commandType='wsclean', processors='max')
    sch.run(check=True)

with w.if_todo('produce_synchrotron_fits'):

    ds = yt.load(pathdir)

    ds.add_field(("gas", "synchrotron_emissivity"), function=_synchrotron_emissivity, units="erg/s/cm**2", sampling_type="local", force_override=True)
    # mylist = ds.derived_field_list
    # for field in mylist:
    #     print(field)
    # sys.exit()

    #width = (50, "kpc")
    resolution = (imsize, imsize)

    # ad = ds.all_data()  # Seleziona l'intero volume del dataset
    # density_data = ad["gas", "synchrotron_emissivity"]
    # print(density_data)
    # sys.exit()

    slice_plot = yt.SlicePlot(ds, "z", ("gas", "synchrotron_emissivity"))
    slice_plot.set_buff_size(resolution)

    frb = slice_plot.frb
    emissivity_data = frb["gas", "synchrotron_emissivity"].v

    hdu = fits.PrimaryHDU(emissivity_data)
    hdu.header['BUNIT'] = 'ERG/S/CM**2'

    hdu.writeto('images/synchrotron_emissivity.fits', overwrite=True)

with w.if_todo('inject_data'):

    emissivity = 'images/synchrotron_emissivity.fits'
    emissivity_hdu = fits.open(emissivity)
    emissivity_data = emissivity_hdu[0].data

    # Convert to Jy
    I_Jy_per_beam = (emissivity_data / 1e-23)

    emissivity_data_reshaped = emissivity_data[np.newaxis, np.newaxis, :, :]

    empty_models = 'images/empty-0*-model.fits'
    input_files = glob.glob(empty_models)

    for model in input_files:
        model_hdu = fits.open(model)
        model_data = model_hdu[0].data

        if emissivity_data_reshaped.shape != model_data.shape:
            raise ValueError(f"Dimensions are not matching for {model}.")

        bmaj = model_hdu[0].header['BMAJ']
        bmin = model_hdu[0].header['BMIN']
        freq = model_hdu[0].header['CRVAL3']

        # beam_area = (np.pi * bmaj * bmin) / (4 * np.log(2))
        # beam_area_pix = np.pi * (bmaj / 2) * (bmin / 2)
        # emissivity_jy_per_beam = emissivity_data_reshaped * beam_area_pix

        model_data[:] = I_Jy_per_beam * freq ** (-(2.5-1)/2) # I add here the frequency dependency

        file_number = os.path.basename(model).split('-')[1]
        output_file = f'models/injected-{file_number}-model.fits'

        hdu = fits.PrimaryHDU(data=model_data, header=model_hdu[0].header)
        hdu.header['BUNIT'] = 'JY/PIXEL'
        hdu.writeto(output_file, overwrite=True)

        model_hdu.close()
        logger.info(f"Injection done in {output_file}.")

    emissivity_hdu.close()

with w.if_todo('predict'):

    logger.info(f'Predicting...')
    sch.add(f'wsclean -predict -name models/injected -channels-out {chout} data/{name}*.MS', log='predict.log', commandType='wsclean', processors='max')
    sch.run(check=True)

    MSs_empty.addcol('CLEAN_DATA', 'MODEL_DATA', log='$nameMS_addcol.log')

    logger.info(f'Producing image with no corruptions...')
    sch.add(f'wsclean -size {imsize} {imsize} -name images/clean -scale {scale} -data-column CLEAN_DATA -weight briggs -1 -circular-beam -niter 100000 -no-update-model-required -mgain 0.6'
        f' -baseline-averaging 10 -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
        log='wsclean-clean.log', commandType='wsclean', processors='max')
    sch.run(check=True)

if not nocorrupt:

    os.system(f'cp {synth}/parsets/empty.skymodel .')

    if recorrupt:
        MSs_empty.deletecol('CORRUPTED_DATA')

        with open('synthetics.walker', 'r') as file:
            lines = file.readlines()

        filtered_lines = [line for line in lines if
                          not (line.startswith('corrupt') or line.startswith('corrupted_image'))]

        with open('synthetics.walker', 'w') as file:
            file.writelines(filtered_lines)

    MSs_empty.addcol('CORRUPTED_DATA', 'MODEL_DATA', log='$nameMS_addcol.log')

    # I need to update the msin of the parset file at each run, in losito there is no way to give it in the command line...
    with open(parset, "r") as file:
        parset_content = file.readlines()

    msin_path = f"data/{name}*.MS"
    new_parset_content = []
    for line in parset_content:
        if line.startswith("msin"):
            new_parset_content.append(f"msin = {msin_path}\n")
        else:
            new_parset_content.append(line)

    with open(parset, "w") as file:
        file.writelines(new_parset_content)

    with w.if_todo('corrupt'):
        logger.info(f'Corrupting visibilities...')
        sch.add(f'losito {parset}')
        sch.run(check=True)

    with w.if_todo('corrupted_image'):
        logger.info(f'Producing image with {corrtype} corruptions...')
        sch.add(f'wsclean -size {imsize} {imsize} -name images/corrupted_{corrtype} -scale {scale} -data-column CORRUPTED_DATA -weight briggs -1 -circular-beam -niter 100000 -no-update-model-required -mgain 0.6'
            f' -baseline-averaging 10 -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
            log='wsclean-corrupted.log', commandType='wsclean', processors='max')
        sch.run(check=True)


