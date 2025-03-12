from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
# import BASE_DIR from settings
from django.conf import settings

import os
def show_exc(exception):
    # return the exception as a string, with file and line number
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    filename_rel = os.path.relpath(filename, os.path.dirname(__file__))
    app_folder = os.path.basename(os.path.dirname(__file__))
    return f'EXCEPTION IN ({filename_rel}, LINE {lineno}): {exception} (APP: {app_folder})'


# Create your views here.
import json, sys
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
warnings.simplefilter("ignore", category=UserWarning)


from astropy.wcs import WCS

from karabo.imaging.imager_base import DirtyImagerConfig
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.imaging.imager_rascil import (
    RascilDirtyImager,
    RascilDirtyImagerConfig,
    RascilImageCleaner,
    RascilImageCleanerConfig,
)
from karabo.imaging.imager_wsclean import (
    WscleanDirtyImager,
    WscleanImageCleaner,
    WscleanImageCleanerConfig,
    create_image_custom_command,
    )
from karabo.simulation.visibility import *
from astropy.coordinates import SkyCoord  # High-level coordinates
from karabo.simulator_backend import SimulatorBackend

from datetime import datetime

import numpy as np

from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.util.plotting_util import get_slices
from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType
from karabo.simulation.telescope import RASCILTelescopes



import time
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz
import astropy.units as u

import uuid

def printlog(msg, t0=None):
    # t0 is the initial timestamp
    if (t0):
        date_for_human = datetime.fromtimestamp(t0).strftime('%Y-%m-%d %H:%M:%S')
        print (f"[{date_for_human} + {time.time()-t0:.2f}sec] {msg}")
    else:
        print (f"[{datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

def optimal_manhattan_point(points):
    # The coordinates x (rows) and y (columns) are separated.
    xs = sorted(p[0] for p in points)
    ys = sorted(p[1] for p in points)
    n = len(points)
    
    # If the number of points is odd, the median is the central value.
    # If it is even, any value between the two central values
    # here we use the average (rounded to integer if the grid is discrete).
    median_x = xs[n // 2] if n % 2 == 1 else (xs[n // 2 - 1] + xs[n // 2]) // 2
    median_y = ys[n // 2] if n % 2 == 1 else (ys[n // 2 - 1] + ys[n // 2]) // 2
    return (median_x, median_y)



def calculate_lst(longitude, latitude, observation_time_local):
    """
    Calculates the Local Sidereal Time (LST) and the equatorial coordinates of the zenith
    for a given location and local time.
    
    Parameters:
        longitude (float): Longitude in degrees (positive for East, negative for West).
        latitude (float): Latitude in degrees (positive for North, negative for South).
        observation_time_local (str): Date and time in 'YYYY-MM-DD HH:MM:SS' format in local time.
    
    Returns:
        dict: Contains the LST in hours, the RA of the zenith in degrees, and the Dec of the zenith in degrees.
    """
    location = EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=300*u.m)
    
    # Convert local time to UTC
    observation_time = Time(observation_time_local, format='iso', scale='utc') - (longitude / 15.0) * u.hour
    
    lst = observation_time.sidereal_time('apparent', longitude*u.deg)
    
    return {
        "LST": lst.to_string(unit=u.hour, sep=':'),
        "RA_zenith": lst.degree,
        "Dec_zenith": latitude
    }


def index(request):
    try:
        telescopes_rascil = get_args(RASCILTelescopes)
        telescopes_oskar = get_args(OSKARTelescopesWithoutVersionType) # + get_args(OSKARTelescopesWithVersionType)
        backend_names = ['OSKAR', 'RASCIL']

        return render(request, 'web/index.html', {'telescopes_oskar': telescopes_oskar, 'telescopes_rascil': telescopes_rascil, 'backend_names': backend_names})
    except Exception as e:
        return render(request, 'exception.html', {'msg': show_exc(e)})
    
def simulation(request):
    error = True
    image1 = None
    image2 = None
    image3 = None

    try:
        t0 = time.time()
        # Get the parameters from the form

        telescope_name = request.POST.get('telescope')
        backend_name = request.POST.get('backend')
        fov = float(request.POST.get('fov', 20))
        N_srcs = int(request.POST.get('N_srcs', 5))
        rms = float(request.POST.get('rms', 0))
        pixel = float(request.POST.get('pixel', 0.5))
        # cfg is a JSON File with the configuration (FILE input type)
        cfg = request.FILES.get('cfg', None)
        if cfg:
            cfg = cfg.read().decode('utf-8')
        else:
            pass

        fout = request.POST.get('fout', f'sim{time.time():.0f}.png')
        if not fout:
            fout = f'sim{time.time():.0f}.png'
        if not fout.endswith('.png'):
            fout += '.png'
    
        maxflux = float(request.POST.get('maxflux', 10))
        observation_time_local = request.POST.get('observation_time_local', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        observation_time_local = datetime.strptime(observation_time_local, '%Y-%m-%d %H:%M:%S')
        cleaned = request.POST.get('cleaned', 'off')
        printlog(cleaned, t0)
        cleaned = cleaned == 'on'

        # add current path to the output file
        # fout = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'simulations' , fout)
        uuid_simulation = str(uuid.uuid4())
        folder_name = os.path.join(settings.BASE_DIR, 'static', 'simulations', uuid_simulation)
        os.makedirs(folder_name, exist_ok=True)
        fout = os.path.join(folder_name, fout)

        if backend_name.lower() == "rascil":
            backend = SimulatorBackend.RASCIL
            telescope_types = get_args(RASCILTelescopes)
            if telescope_name in telescope_types:
                telescope = Telescope.constructor(telescope_name, backend=backend)
            else:
                printlog (f"Telescope {telescope_name} not found. Loading LOFAR as default. Values accepted: {telescope_types}", t0)
                telescope = Telescope.constructor('LOFAR', backend=backend)

        else:
            backend = SimulatorBackend.OSKAR
            telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
            if (telescope_name in telescope_types):
                telescope = Telescope.constructor(telescope_name, backend=backend)
            else:
                printlog (f"Telescope {telescope_name} not found. Loading LOFAR as default. Values accepted: {telescope_types}", t0)
                telescope = Telescope.constructor("LOFAR", backend=backend)
        printlog (f"Telescope loaded: {telescope.name.upper()}", t0)

        printlog (f"Generating {N_srcs} random sources", t0)
        limit = fov/2
        result = calculate_lst(telescope.centre_longitude, telescope.centre_latitude, observation_time_local.strftime("%Y-%m-%d %H:%M:%S"))
        x0 = result["RA_zenith"]
        y0 = result["Dec_zenith"]
        # Random sources [x0, y0, flux]; x0, y0 is the limited between -20, 20. Flux is limited between [0.1, 10] 
        sky_data = np.random.rand(N_srcs, 3) * 2*limit - limit
        sky_data[:, 2] = ((sky_data[:, 2] + limit) / (2*limit)) * (maxflux - 0.1) + 0.1
        sky_data[:, 0] += x0
        sky_data[:, 1] += y0

        observation = Observation(
            start_frequency_hz=1e6,
            start_date_and_time=observation_time_local,
            phase_centre_ra_deg = x0,
            phase_centre_dec_deg = y0,
        )

        imaging_cellsize = (pixel *  u.arcsec).to(u.deg).value
        imaging_npixel = 2048
        maxflux = np.max(sky_data[:, 2])

        print (f"Imaging cellsize: {imaging_cellsize} deg")
        print (f"Imaging npixel: {imaging_npixel}")


        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ['RA---SIN', 'DEC--SIN']
        wcs.wcs.crval = [x0, y0]
        wcs.wcs.crpix = [imaging_npixel//2, imaging_npixel//2]
        wcs.wcs.cdelt = [-imaging_cellsize, +imaging_cellsize]
        wcs.wcs.radesys = 'ICRS'
        wcs.wcs.equinox = 2000.0
        
        sky = SkyModel(wcs=wcs)
        sky.add_point_sources(sky_data)
        sky.explore_sky([x0, y0],wcs=wcs, filename=fout.replace(".png", "_sources.png"), xlabel='RA', ylabel='DEC', vmin=0)

        # run a single simulation with the provided configuration 
        simulation = InterferometerSimulation()

        # Get current path
        current_path = os.path.dirname(os.path.realpath(__file__))
        visibility_path = os.path.join(current_path, 'MSDIRS', f"{uuid_simulation}.MS") # path to the visibility file; if you use the WSCLEAN backend, the path must be absolute
        simulation.run_simulation(telescope, sky, observation, visibility_path=visibility_path)
        visibilities = Visibility(visibility_path)


        # Create the dirty image
        imaging_cellsize = (pixel *  u.arcsec).to(u.deg).value
        imaging_npixel = 2048

        if backend_name.lower() == "rascil":
            config = RascilDirtyImagerConfig(
                    imaging_npixel=imaging_npixel,
                    imaging_cellsize=imaging_cellsize,
                )
            imager = RascilDirtyImager(config)
            
        else:
            config = OskarDirtyImagerConfig(
                    imaging_phase_centre=SkyCoord(ra=x0, dec=y0, unit="deg"),
                    imaging_npixel=imaging_npixel,
                    imaging_cellsize=imaging_cellsize,
                )
            imager = OskarDirtyImager(config)

        dirty = imager.create_dirty_image(visibilities)
        if rms > 0:
            printlog (f"Adding noise to the image with RMS: {rms}", t0)
            noise_matrix = np.random.normal(0, rms, dirty.data[0].shape)
            dirty.data[0] += noise_matrix


        printlog(f"Saving dirty image to {fout}", t0)
        dirty.plot(title=f"Dirty image {backend_name.upper()} ({telescope.name.upper()})", filename=fout, wcs_enabled=True, xlabel='RA', ylabel='DEC', vmax=maxflux * 1.05, vmin=0)

        printlog (f"Sky Model with {len(sky_data)} sources", t0)
        printlog (f"Optimal phase center: {x0}, {y0}", t0)

        sky_data = np.array(sky_data)

        # Convert sky_data to json
        sky_data = sky_data.tolist()
        sources = []
        for i, src in enumerate(sky_data):
            sources.append({"ra": src[0], "dec": src[1], "flux": src[2], "mute":False, "name":f"source_{i:03}"})
        json_data = {"sources": sources, "phase_center": {"ra": x0, "dec": y0}, "observation_date": observation_time_local.strftime("%Y-%m-%d %H:%M:%S")}
        # if (cfg is None):
        #     with open(fout.replace(".png", "_sources.json"), "w") as f:
        #         f.write(json.dumps(json_data, indent=4))
        #     printlog (f"Saved sources to {fout.replace('.png', '_sources.json')}", t0)
        printlog(cleaned, t0)
        if cleaned:
            printlog ("Cleaning the image...", t0)
            if backend_name.lower() in ["rascil", "all"]:
                try:
                    deconvolved, restored, residual = RascilImageCleaner(
                        RascilImageCleanerConfig(
                            imaging_npixel=imaging_npixel,
                            imaging_cellsize=imaging_cellsize,
                            ingest_vis_nchan=1,
                            # clean_nmajor=1,
                            # clean_algorithm="mmclean",
                            # clean_scales=[10, 30, 60],
                            clean_threshold=0.12e-3,
                            clean_nmoment=5,
                            # clean_psf_support=640,
                            clean_restored_output="integrated",
                            use_dask=False,
                        )
                    ).create_cleaned_image(visibilities)
                    restored.plot(title=f"Cleaned image (RASCIL) {backend_name.upper()} ({telescope.name.upper()})", filename=fout.replace(".png", "_cleaned.png"), wcs_enabled=True, xlabel='RA', ylabel='DEC')
                    printlog (f"Saved cleaned image to {fout.replace('.png', '_cleaned.png')}", t0)
                except Exception as e:
                    printlog (f"Error cleaning the image: {e}", t0)
                    printlog ("Using WSCLEAN instead", t0)
                    config = WscleanImageCleanerConfig(imaging_npixel=imaging_npixel, imaging_cellsize=imaging_cellsize)
                    cleaner = WscleanImageCleaner(config)
                    cleaned = cleaner.create_cleaned_image(visibilities, output_fits_path=fout.replace(".png", "_wccleaned.fits"))
                    cleaned.plot(title=f"Cleaned image (WSCLEAN) {args.backend.upper()} ({telescope.name.upper()})", filename=fout.replace(".png", "_cleaned.png"), wcs_enabled=True, xlabel='RA', ylabel='DEC')
                    printlog (f"Saved cleaned image to {fout.replace('.png', '_cleaned.png')}", t0)
            if backend_name.lower() in ["oskar", "all", "wsclean"]:
                printlog ("Cleaning not supported for OSKAR, using WSCLEAN", t0)
                config = WscleanImageCleanerConfig(imaging_npixel=imaging_npixel, imaging_cellsize=imaging_cellsize)
                cleaner = WscleanImageCleaner(config)
                cleaned = cleaner.create_cleaned_image(visibilities, output_fits_path=fout.replace(".png", "_wccleaned.fits"))
                cleaned.plot(title=f"Cleaned image (WSCLEAN) {backend_name.upper()} ({telescope.name.upper()})", filename=fout.replace(".png", "_cleaned.png"), wcs_enabled=True, xlabel='RA', ylabel='DEC')
                printlog (f"Saved cleaned image to {fout.replace('.png', '_cleaned.png')}", t0)

        image1 = os.path.join('static','simulations', uuid_simulation, os.path.basename(fout))
        image2 = os.path.join('static','simulations', uuid_simulation, os.path.basename(fout.replace(".png", "_sources.png")))
        image3 = os.path.join('static','simulations', uuid_simulation, os.path.basename(fout.replace(".png", "_cleaned.png")))


        return JsonResponse({'error': False, 'error_msg': 'Ok', 'image1': image1, 'image2': image2, 'image3': image3})
    except Exception as e:
        return JsonResponse({'error': True, 'error_msg': show_exc(e)})


    return JsonResponse({'error': True, 'error_msg': 'Ok', 'image1': image1, 'image2': image2, 'image3': image3})
    return HttpResponse("Hello, world. You're at the simulation page.")
    return render(request, 'web/simulation.html')