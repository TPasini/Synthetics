#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2024 - Thomas Pasini
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys, os, glob, argparse
import numpy as np
from synth_libs import lib_ms, lib_log, lib_util
import warnings
import astropy
from astropy.cosmology import FlatLambdaCDM
import yt
from yt.units import gauss, g, cm, Hz, W, m, Jy, erg, s, G
from astropy.io import fits
import math
import bdsf
from itertools import combinations
from synth_libs.lib_util import create_region
import h5py

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
scale = '3asec'

##################################################################

#This is just a very general and simple synchrotron function where we assume constant density of relativistic electrons.
def _synchrotron_base(field, data):
    #n_e = data["gas", "density"]/9.11e-28 # Converting density into numeric density #TODO this is gas density, we need the relativistic electron density
    n_e = 1e-5 * cm**-3
    p = 2.5  # Power-law index (2.5 is quite typical...)
    B_perp = abs(data["gas", "magnetic_field_y"]*1e-6) #We need this 1e-6 factor because for some reason it's missing when read from within the function...
    emissivity = n_e * (B_perp** ((p+1)/2)) / G**(7/4) * erg * s**-1 * cm**-2 * cm**3 #We add the frequency dependency later because it's not in the simulation
    return emissivity

def _synchrotron_emissivity(field, data):
    field_list = data.ds.field_list

    if ("gas", "density") in field_list:
        dens_data = data["gas", "density"]
    else:
        dens_key = next((f for f in field_list if f[0] == "stream" and "dens" in f[1].lower()), None)
        if dens_key is None: raise ValueError("Campo densità non trovato!")
        dens_data = data[dens_key]

    if ("gas", "magnetic_field_y") in field_list:
        by_data = data["gas", "magnetic_field_y"]
    else:
        by_key = next((f for f in field_list if f[0] == "stream" and f[1] in ["By", "magnetic_field_y", "B_y"]), None)
        if by_key is None: raise ValueError("Campo magnetico Y non trovato!")
        by_data = data[by_key]

    n_e = (dens_data / 9.1e-28) * 1e-4
    p = 2.5
    B_perp = abs(by_data * 1e-6)
    emissivity = n_e * B_perp ** 2 / G ** 2 * erg * s ** -1 * cm ** -2 * cm ** 3 * g ** -1
    return emissivity

#The following function is written in a way that relevant fields can be read from a more generic simulation.
def _synchrotron_emissivity(field, data):
    z = 0.1
    cdd = 2.82e-30 * (1 + z) ** 3  # Density factor -> g/cm^3
    cv = 2.51e9
    cb = np.sqrt(cdd * 4 * np.pi) * cv * (1 + z) ** 0.5  # Magnetic field factor -> Gauss

    n_e_rel = 1e-6 * (cm ** -3)

    field_list = data.ds.field_list

    if ("gas", "density") in field_list:
        dens_raw = data["gas", "density"]
    else:
        dens_key = next((f for f in field_list if f[0] == "stream" and "dens" in f[1].lower()), None)
        if not dens_key: raise ValueError("Campo Density non trovato!")
        dens_raw = data[dens_key]

    b_components = {}
    for ax in ['x', 'y', 'z']:
        standard_f = ("gas", f"magnetic_field_{ax}")
        if standard_f in field_list:
            b_components[ax] = data[standard_f]
        else:
            stream_f = next((f for f in field_list if f[0] == "stream" and f[1].lower() == f"b{ax}"), None)
            if not stream_f: raise ValueError(f"Campo B{ax} non trovato!")
            b_components[ax] = data[stream_f]

    rho_g_cm3 = dens_raw * cdd

    Bx_G = b_components['x'] * cb
    By_G = b_components['y'] * cb
    Bz_G = b_components['z'] * cb

    B2 = (Bx_G ** 2 + By_G ** 2 + Bz_G ** 2)

    emissivity = n_e_rel * B2 * erg * s ** -1 * cm

    return emissivity

def convert_coordinates(lat_deg, lon_deg):
    lat_rad = math.radians(lat_deg)
    lon_rad = math.radians(lon_deg)
    return lat_rad, lon_rad

##################################################################

parser = argparse.ArgumentParser(description='Create synthetic .MS file(s) which contain the mock visibilities of the synchrotron emission derived from a simulation cube, as detected by the LOFAR telescope.')
parser.add_argument('-p', '--path', dest='path', action='store', default='', type=str, help='Path where to look for the simulation cube.')
parser.add_argument('-radec', '--radec', dest='radec', nargs=2, type=float, default=None, help='RA/DEC of the phase centre of the generated .MS file.')
parser.add_argument('--name', dest='name', type=str, default='simulated', help='Name of the generated .MS file.')
parser.add_argument('--begin', dest='begin', type=float, default='61400', help='MJD starting time of the generated .MS file.')
parser.add_argument('--duration', dest='duration', type=float, default=1, help='Duration of the generated .MS file in hours. Default is 1 hour.')
parser.add_argument('--station', dest='station', type=str, default='LBA', help='Whether to use LBA or HBA. Default is LBA.')
parser.add_argument('--version', dest='version', type=int, default=1, help='Whether to use LOFAR or LOFAR 2.0. Default is LOFAR. Use 2 for LOFAR 2.0')
parser.add_argument('--freqrange', dest='freqrange', nargs=2, type=float, default=(42e6, 66e6), help='Frequency range of the generated .MS file in Hz.')
parser.add_argument('--chanpersb', dest='chanpersb', type=int, default=1, help='Number of channels per sub-band. Default is 1. A maximum of 4 is supported by LoSiTo.')
parser.add_argument('--chout', dest='chout', type=int, default=6, help='Number of output channel images in WSClean. Default is 6.')
parser.add_argument('--imsize', dest='imsize', type=int, default=1024, help='Pixel size of images. Default is 1024.')
parser.add_argument('--nocorrupt', dest='nocorrupt', action='store_true', help='Whether to corrupt the dataset with ionospheric delays.')
parser.add_argument('--corruption_type', dest='corrtype', action='store', nargs='+', default="all", type=str, help='Type of corruption to apply to the dataset. Can be set to tec, fr, clock, polmisalign, beamcorrupt, noise or all, separated by spaces. Default is all.')
parser.add_argument('--recorrupt', dest='recorrupt', action='store_true', help='Use this if you just want to change the type of corruption to apply, without re-running everything else.')
parser.add_argument('--skymodel_bdsf', dest='skymodel_bdsf', action='store_true', help='Use PyBDSF to create a sky model for corruptions, instead of predicting from the initial .fits image. Only advanced users.')


args = parser.parse_args()
pathdir = args.path
deg_coords = args.radec
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
skymodel_bdsf = args.skymodel_bdsf

if not os.path.exists(name):
    os.makedirs(name)
os.chdir(name)

logger_obj = lib_log.Logger('synthetics')
logger = lib_log.logger
sch = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('synthetics.walker')
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)

#Since we want to give coordinates in degrees, we need to convert to radians for synthms
coords = convert_coordinates(deg_coords[0], deg_coords[1])

valid_elements = ["tec", "fr", "clock", "polmisalign", "beamcorrupt", "noise", "all"]

valid_combinations = set()
for r in range(1, len(valid_elements) + 1):  # Combina da 1 elemento fino a tutti
    valid_combinations.update(tuple(sorted(combo)) for combo in combinations(valid_elements, r))

# if tuple(sorted(corrtype)) not in valid_combinations:
#     logger.error(f'In {corrtype} there is an unknown corruption type. Possible values are combinations of tec, fr, clock, polmisalign, beamcorrupt, noise, or all.')
#     sys.exit()
# else:
#     logger.info(f'Corruptions for {corrtype} will be applied to the observation.')

if not pathdir:
    logger.error('Provide a path (-p) where to look for the simulation cube.')
    sys.exit()

if not coords:
    logger.error('Provide RA and DEC in degrees of the phase centre of the .MS file you want to generate.')
    sys.exit()

if chanpersb > 4:
    logger.error('LoSiTo supports a maximum of 4 channels per subband.')
    sys.exit()

if ('tec' in corrtype or 'fr' in corrtype) and skymodel_bdsf == False:
    logger.error('Currently only direction-independent corruptions are supported when predicting from a .fits image. Either remove direction-dependent corruptions or specify --skymodel_bdsf in the command line.')
    sys.exit()

with w.if_todo('cleaning'):
    logger.info('Preparing the environment...')
    lib_util.check_rm('data')
    lib_util.check_rm('images')
    lib_util.check_rm('models')
    lib_util.check_rm('skymodels')
    lib_util.check_rm('parsets')
    os.makedirs('data')
    os.makedirs('images')
    os.makedirs('models')
    os.makedirs('skymodels')
    os.makedirs('parsets')

with w.if_todo('generate_MS'):
    logger.info('Generating empty .MS files...')
    sch.add(f'synthms --name data/{name} --tobs {duration} --ra {coords[0]} --dec {coords[1]} --station {station} '
          f'--lofarversion {version} --minfreq {freqrange[0]} --maxfreq {freqrange[1]} --chanpersb {chanpersb}', log='generateMS.log', commandType='general', processors='max')
    sch.run(check=False)

logger.info('Opening MS files...')
MSs_empty = lib_ms.AllMSs(glob.glob(f'data/{name}*.MS'), sch, check_flags=False)

with w.if_todo('empty_image'):
    logger.info('Producing empty image from initial dataset...')
    sch.add(f'wsclean -size {imsize} {imsize} -name images/empty -scale {scale} -data-column DATA -weight briggs -0.3 -circular-beam -niter 100000 -no-update-model-required -mgain 0.6'
          f' -baseline-averaging 10 -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
          log='wsclean-empty.log', commandType='wsclean', processors='max')
    sch.run(check=True)

with w.if_todo('produce_injection_fits'):
    try:
        ds = yt.load(pathdir)
        _ = ds.field_list

    except Exception as e:
        data_dict = {}
        def extract_3d_datasets(name, obj):
            if isinstance(obj, h5py.Dataset):
                if len(obj.shape) == 3:
                    key = name.split('/')[-1]
                    data_dict[key] = obj[:]

        with h5py.File(pathdir, "r") as f:
            f.visititems(extract_3d_datasets)

        if not data_dict:
            logger.error("No 3D arrays found in the input file.")
            sys.exit()

        logger.info(f"Found the following fields: {list(data_dict.keys())}.")

        grid_shape = list(data_dict.values())[0].shape
        bbox = np.array([[0, 1], [0, 1], [0, 1]])

        ds = yt.load_uniform_grid(data_dict, grid_shape, bbox=bbox)

    ds.add_field(("gas", "synchrotron_emissivity"), function=_synchrotron_emissivity, units="erg/s/cm**2",
                 sampling_type="local", force_override=True)

    resolution = (imsize, imsize)

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
    sch.add(f'wsclean -predict -name models/injected -channels-out {chout} data/{name}*.MS', log='predict_raw.log', commandType='wsclean', processors='max')
    sch.run(check=True)

    MSs_empty.addcol('CLEAN_DATA', 'MODEL_DATA', log='$nameMS_addcol.log')

with w.if_todo('clean_image'):
    logger.info(f'Producing image with no corruptions...')
    sch.add(f'wsclean -size {imsize} {imsize} -name images/clean -scale {scale} -data-column CLEAN_DATA -weight briggs -0.3 -circular-beam -niter 100000 -no-update-model-required -mgain 0.6'
        f' -baseline-averaging 4 -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
        log='wsclean-clean.log', commandType='wsclean', processors='max')
    sch.run(check=True)

with w.if_todo('predict_back'):
    logger.info(f'Predict convolved data back into model...')
    sch.add(f'wsclean -predict -name images/clean -channels-out {chout} data/{name}*.MS', log='predict_convolved.log', commandType='wsclean', processors='max')
    sch.run(check=True)

if not nocorrupt:

    if not nocorrupt and not corrtype or corrtype == ['all']:
        corrtype = ['tec', 'fr', 'clock', 'polmisalign', 'beam', 'noise', 'bandpass']

    corr_list = corrtype_str = "_".join(corrtype)

    if recorrupt:
        MSs_empty.deletecol('CORRUPTED_DATA')

        with open('synthetics.walker', 'r') as file:
            lines = file.readlines()

        filtered_lines = [line for line in lines if
                          not (line.startswith('corrupt') or line.startswith('corrupted_image'))]

        with open('synthetics.walker', 'w') as file:
            file.writelines(filtered_lines)

    MSs_empty.addcol('CORRUPTED_DATA', 'MODEL_DATA', log='$nameMS_addcol.log')

    msin_path = f"data/{name}*.MS"

    if skymodel_bdsf:
        with w.if_todo('create_skymodel'):
        # Run PyBDSF to get a good sky model to use for corruptions
            img = bdsf.process_image('images/clean-MFS-image.fits', rms_map=False, mean_map='zero', atrous_do = True)
            img.write_catalog(outfile=f'skymodels/{corr_list}.skymodel', catalog_type='gaul', format='bbs', clobber=True)
            img.export_image(outfile=f'skymodels/{corr_list}.fits', img_type='gaus_model', clobber=True)

    with w.if_todo('create_region'):
        with fits.open('images/clean-MFS-image.fits') as hdul:
            header = hdul[0].header

        naxis = header['NAXIS1']  # Image width in pixels

        ra_center = header['CRVAL1']  # Central RA
        dec_center = header['CRVAL2']  # Central Dec

        extent = abs(header['CDELT1'] * naxis / 2.0) # We need half the total extent since we start from the centre

        region_str = create_region(ra_center, dec_center, extent*np.sqrt(2), extent, shape='polygon')
        patch_str = f'point({ra_center}, {dec_center})  # point=boxcircle text={{Patch1}}'
        region_str += patch_str
        with open('region.reg', 'w') as f:
            f.write(region_str)


    with open(f'parsets/{corr_list}.parset', "w") as file:
        new_parset_content = []
        new_parset_content.append(f"msin = {msin_path}\n")
        if skymodel_bdsf:
            new_parset_content.append(f"skymodel = skymodels/{corr_list}.skymodel\n")
        else:
            new_parset_content.append(f"skymodel = images/clean-MFS-image.fits\n")
            new_parset_content.append(f"regions = region.reg\n")
        new_parset_content.append(f"\n")

        corruption_parset_content = []
        if 'tec' in corrtype:
            corruption_parset_content.append(f"[tec]\n")
            corruption_parset_content.append(f"operation = TEC\n")
            corruption_parset_content.append(f"method = turbulence\n")
            corruption_parset_content.append(f"\n")
        if 'fr' in corrtype:
            corruption_parset_content.append(f"[faraday]\n")
            corruption_parset_content.append(f"operation = FARADAY\n")
            corruption_parset_content.append(f"\n")
        if 'clock' in corrtype:
            corruption_parset_content.append(f"[clock]\n")
            corruption_parset_content.append(f"operation = CLOCK\n")
            corruption_parset_content.append(f"\n")
        if 'polmisalign' in corrtype:
            corruption_parset_content.append(f"[polmisalign]\n")
            corruption_parset_content.append(f"operation = POLMISALIGN\n")
            corruption_parset_content.append(f"\n")
        if 'beam' in corrtype:
            corruption_parset_content.append(f"[beam]\n")
            corruption_parset_content.append(f"operation = BEAM\n")
            corruption_parset_content.append(f"mode = default\n")
            corruption_parset_content.append(f"\n")

        corruption_parset_content.append(f"[predict]\n")
        corruption_parset_content.append(f"operation = PREDICT\n")
        corruption_parset_content.append(f"outputcolumn = CORRUPTED_DATA\n")
        corruption_parset_content.append(f"resetWeights = True\n")
        if skymodel_bdsf:
            corruption_parset_content.append(f"predictType = h5parmpredict\n")
        else:
            corruption_parset_content.append(f"predictType = idgpredict\n")
        corruption_parset_content.append(f"\n")

        if 'noise' in corrtype:
            corruption_parset_content.append(f"[noise]\n")
            corruption_parset_content.append(f"operation = NOISE\n")
            corruption_parset_content.append(f"outputcolumn = CORRUPTED_DATA\n")
            corruption_parset_content.append(f"\n")
        if 'bandpass' in corrtype:
            corruption_parset_content.append(f"[bandpass]\n")
            corruption_parset_content.append(f"operation = BANDPASS\n")
            corruption_parset_content.append(f"method = ms\n")
            corruption_parset_content.append(f"\n")

        file.writelines(new_parset_content)
        file.writelines(corruption_parset_content)

    with w.if_todo('corrupt'):
        logger.info(f'Corrupting visibilities...')
        sch.add(f'losito parsets/{corr_list}.parset')
        sch.run(check=True)

    with w.if_todo('corrupted_image'):
        logger.info(f'Producing image with {corrtype} corruptions...')
        sch.add(f'wsclean -size {imsize} {imsize} -name images/corrupted_{corr_list} -scale {scale} -data-column CORRUPTED_DATA -weight briggs -0.3 -circular-beam -niter 100000 -mgain 0.6'
            f' -join-channels -fit-spectral-pol 3 -channels-out {chout} data/{name}*.MS',
            log='wsclean-corrupted.log', commandType='wsclean', processors='max')
        sch.run(check=True)

logger.info('Done.')