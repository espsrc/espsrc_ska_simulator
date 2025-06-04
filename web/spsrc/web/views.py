from django.http import HttpResponse, JsonResponse
from django.shortcuts import redirect, render
from django.conf import settings
from django.views.decorators.csrf import csrf_exempt
from matplotlib.colors import PowerNorm
from datetime import timedelta
import astropy.coordinates as acoord


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


class Source:
    def __init__(self, ra, dec, I, Q=0 * u.Jy, U=0 * u.Jy, V=0 * u.Jy, ref_freq=0 * u.Hz, spec_index=0, rot_meas=0 * u.rad/(u.m**2), 
                 major_axis = 0*u.arcsec, minor_axis = 0*u.arcsec, pa=0*u.arcsec,  true_redshift=0, obs_redshift=0, obj_id = None):
        # Initialize the source with its parameters, checking units

        list_of_units = [u.deg, u.deg, u.Jy, u.Jy, u.Jy, u.Jy, u.Hz, u.rad/(u.m**2), u.arcsec, u.arcsec, u.arcsec]
        list_of_values = [ra, dec, I, Q, U, V, ref_freq, rot_meas, major_axis, minor_axis, pa]
        for i, unit in enumerate(list_of_units):
            if not isinstance(list_of_values[i], u.Quantity):
                list_of_values[i] = list_of_values[i] * unit
                print (f"Adding unit {unit} to value {list_of_values[i]}")
            if list_of_values[i].unit != unit:
                try:
                    list_of_values[i] = list_of_values[i].to(unit)
                except u.UnitConversionError:
                    raise ValueError(f"Value {list_of_values[i]} does not have the correct unit {unit}")
        self.ra = list_of_values[0]
        self.dec = list_of_values[1]
        self.I = list_of_values[2]
        self.Q = list_of_values[3]
        self.U = list_of_values[4]
        self.V = list_of_values[5]
        self.ref_freq = list_of_values[6]
        self.spec_index = spec_index
        self.rot_meas = list_of_values[7]
        self.major_axis = list_of_values[8]
        self.minor_axis = list_of_values[9]
        self.pa = list_of_values[10]
        self.true_redshift = true_redshift
        self.obs_redshift = obs_redshift
        self.obj_id = obj_id

    @staticmethod
    def from_name(name):
        source = acoord.get_icrs_coordinates(name)
        if source is None:
            raise ValueError(f"Source {name} not found")
        return Source(source.ra, source.dec, 1 * u.Jy)
    
    def to_json(self):
        # Convert the source to a JSON serializable dictionary
        return {
            "ra": self.ra.to(u.deg).value,
            "dec": self.dec.to(u.deg).value,
            "I": self.I.to(u.Jy).value,
            "Q": self.Q.to(u.Jy).value,
            "U": self.U.to(u.Jy).value,
            "V": self.V.to(u.Jy).value,
            "ref_freq": self.ref_freq.to(u.Hz).value,
            "spec_index": self.spec_index,
            "rot_meas": self.rot_meas.value,
            "major_axis": self.major_axis.to(u.arcsec).value,
            "minor_axis": self.minor_axis.to(u.arcsec).value,
            "pa": self.pa.to(u.arcsec).value,
            "true_redshift": self.true_redshift,
            "obs_redshift": self.obs_redshift,
        }
    
    def __str__(self):
        # Return a string representation of the source (only non-zero values)
        str2print = f"Source(ra={self.ra}, dec={self.dec}, I={self.I}"

        json_values = self.to_json()
        for key, value in json_values.items():
            if key not in ['ra', 'dec', 'I'] and value != 0:
                str2print += f", {key}={value}"
        str2print += ")"
        return str2print
    
    def to_sky_model(self, reduced_form=False):
        # Convert the source to a SkyModel object
        if reduced_form:
            return (self.ra.value, self.dec.value, self.I.value)
        else:
            return( self.ra.value, self.dec.value, self.I.value, 
                self.Q.value, self.U.value, self.V.value,
                self.ref_freq.value, self.spec_index,
                self.rot_meas.value, self.major_axis.value,
                self.minor_axis.value, self.pa.value,
                self.true_redshift, self.obs_redshift)
        
    def coords(self, frame='icrs'):
        # Return the coordinates of the source
        return SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame=frame)

    def get_best_observation_time(self, telescope:Telescope, date=None):
        """
        Returns the local time at which an object with a given RA/Dec culminates (best observation time).
        
        Parameters:
        - ra_hours: Right Ascension in hours (float)
        - dec_degrees: Declination in degrees (float)
        - lat_deg: Observer's latitude in degrees (float)
        - lon_deg: Observer's longitude in degrees (float, positive to the East)
        - elevation_m: Altitude above sea level (optional)
        - date: Date as a string 'YYYY-MM-DD' (optional, defaults to today if not provided)
        - timezone_offset: Time difference relative to UTC (e.g., -6 for CDMX)

        Returns:
        - Best time.
        """

        if date is None:
            date = datetime.now().strftime('%Y-%m-%d')

        coord = SkyCoord(ra=self.ra, dec=self.dec)
        location = EarthLocation(lat=telescope.centre_latitude * u.deg, lon=telescope.centre_longitude * u.deg, height=telescope.centre_altitude * u.m)

        midnight = Time(f"{date} 00:00:00") + 12*u.hour  # Mediodía UTC
        delta = timedelta(minutes=1)

        best_time = None
        max_alt = -90

        for minutes in range(-360, 360):
            current_time = midnight + minutes * u.minute
            altaz = coord.transform_to(AltAz(obstime=current_time, location=location))
            if altaz.alt.deg > max_alt:
                max_alt = altaz.alt.deg
                best_time = current_time

        return best_time


def printlog(fname, *args):
    # Print to console and log file
    print (f'[{datetime.now()}]',*args)
    with open(fname, 'a') as f:
        print(f'[{datetime.now()}]', *args, file=f)
        f.flush()
        os.fsync(f.fileno())
        f.close()

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
    
    # Convert local time to UTC
    observation_time = Time(observation_time_local, format='iso', scale='utc') - (longitude / 15.0) * u.hour
    
    lst = observation_time.sidereal_time('apparent', longitude*u.deg)
    
    return {
        "LST": lst.to_string(unit=u.hour, sep=':'),
        "RA_zenith": lst.degree,
        "Dec_zenith": latitude
    }


def index(request, json_data=None):
    try:
        telescopes_rascil = get_args(RASCILTelescopes)
        telescopes_oskar = get_args(OSKARTelescopesWithoutVersionType) # + get_args(OSKARTelescopesWithVersionType)
        backend_names = ['OSKAR', 'RASCIL']

        return render(request, 'web/index.html', {'telescopes_oskar': telescopes_oskar, 'telescopes_rascil': telescopes_rascil, 'backend_names': backend_names})
    except Exception as e:
        return render(request, 'exception.html', {'msg': show_exc(e)})

@csrf_exempt
def simulation(request):
    return simulation_background(request)

@csrf_exempt
def simulation_background(request, uuid_simulation=None):
    error = True
    image1 = None
    image2 = None
    image3 = None

    try:
        t0 = time.time()
        # Get the parameters from the form
        telescope_name = request.POST.get('telescope')
        backend_name = request.POST.get('backend')
        fov = float(request.POST.get('fov', 0))
        ra_list = request.POST.getlist('ra[]', None)
        dec_list = request.POST.getlist('dec[]', None)
        flux_list = request.POST.getlist('flux[]', None)
        pos_list = request.POST.getlist('position[]', None)
        frequency = float(request.POST.get('freq', 1)) * u.MHz
        delta_freq = float(request.POST.get('delta_freq', 1)) * u.MHz
        n_channels = int(request.POST.get('channels', 4))
        seconds  = int(request.POST.get('exptime', 60))



        if fov == 0:
            # Convert frequency to wavelength
            wavelength = frequency.to(u.m, equivalencies=u.spectral())
            fov = (1.22 * wavelength / (13 * u.m)) * u.rad
        else:
            fov = fov * u.deg


        



        rms = float(request.POST.get('rms', 0))
        pixel_size = float(request.POST.get('pixel', 0.5))

        pixel_size = pixel_size * u.arcsec
        tofits = (request.POST.get('tofits', 'off') == 'on')
        tofits = True

    
        maxflux = float(request.POST.get('maxflux', 10))
        observation_time_local = request.POST.get('observation_time_local', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        observation_time_local = datetime.strptime(observation_time_local, '%Y-%m-%d %H:%M:%S')
        cleaned = request.POST.get('cleaned', 'off')
        cleaned = cleaned == 'on'

        if uuid_simulation is None:
            uuid_simulation = str(uuid.uuid4())

        log_file = os.path.join(settings.BASE_DIR, 'static', 'simulations', uuid_simulation, f'{uuid_simulation}.log')
        folder_name = os.path.join(settings.BASE_DIR, 'static', 'simulations', uuid_simulation)
        dirty_path = os.path.join(folder_name, "dirty.png")
        sources_path = os.path.join(folder_name, "sources.png")
        cleaned_path = os.path.join(folder_name, "cleaned.png")
        cfg_path = os.path.join(folder_name, "cfg.json")
        os.makedirs(folder_name, exist_ok=True)

        if backend_name.lower() == "rascil":
            backend = SimulatorBackend.RASCIL
            telescope_types = get_args(RASCILTelescopes)
            if telescope_name in telescope_types:
                telescope = Telescope.constructor(telescope_name, backend=backend)
            else:
                printlog (log_file, f"Telescope {telescope_name} not found. Loading LOFAR as default. Values accepted: {telescope_types}")
                telescope = Telescope.constructor('LOFAR', backend=backend)

        else:
            backend = SimulatorBackend.OSKAR
            telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
            if (telescope_name in telescope_types):
                telescope = Telescope.constructor(telescope_name, backend=backend)
            else:
                printlog (log_file, f"Telescope {telescope_name} not found. Loading LOFAR as default. Values accepted: {telescope_types}")
                telescope = Telescope.constructor("LOFAR", backend=backend)
        printlog (log_file, f"Telescope loaded: {telescope.name.upper()}")

        limit = pixel_size.to(u.deg).value * 2048
        result = calculate_lst(telescope.centre_longitude, telescope.centre_latitude, observation_time_local.strftime("%Y-%m-%d %H:%M:%S"))
        x0 = result["RA_zenith"]
        y0 = result["Dec_zenith"]
        # Random sources [x0, y0, flux]; x0, y0 is the limited between -20, 20. Flux is limited between [0.1, 10] 
        if not (ra_list and dec_list and flux_list):
            N_srcs = int(request.POST.get('N_srcs', 5))
            if (N_srcs == 1):
                limit = 0.001
            printlog (log_file, f"Generating {N_srcs} random sources")
            ra_list = np.random.uniform(x0 - limit, x0 + limit, N_srcs)
            dec_list = np.random.uniform(y0 - limit, y0 + limit, N_srcs)
            flux_list = np.random.uniform(0.1, maxflux, N_srcs)
            pos_list = np.random.randint(1, 10, N_srcs)

        else:
            printlog (log_file, f"Loading {N_srcs} sources from the form")
            N_srcs = len(ra_list)
            ra_list = [float(x) for x in ra_list]
            dec_list = [float(x) for x in dec_list]
            flux_list = [float(x) for x in flux_list]
            pos_list = [int(x) for x in pos_list]

        sources = []
        sources_list = []
        for idx in range(len(ra_list)):
            source = Source(ra=ra_list[idx], 
                                    dec=dec_list[idx],
                                    I=flux_list[idx], 
                                    obj_id=idx, spec_index=0, 
                                    ref_freq=frequency, 
                                    true_redshift=0, 
                                    obs_redshift=0)
            sources.append(source)
            sources_list.append(source.to_sky_model(reduced_form=True))

        try:
            printlog (log_file, f"Creating SkyModel with {len(sources_list)} sources")
            sources_list = np.array(sources_list)
            sky = SkyModel()
            sky.add_point_sources(sources_list)
            printlog (log_file, f"SkyModel created with {len(sky.sources)} sources")
            sky.explore_sky([x0, y0], filename=sources_path, vmin=0, vmax=1.05 * np.max(sources_list[:, 2]),cfun=np.abs)
        except Exception as e:

            printlog (log_file, f"Error creating SkyModel: {show_exc(e)}")
            return JsonResponse({'error': True, 'error_msg': show_exc(e)})
        
        printlog (log_file, f"Sky Model with {len(sources)} sources")
        printlog (log_file, f"Optimal phase center: {x0}, {y0}")
        
        start_freq = frequency - (n_channels/2) * delta_freq
        # imaging_cellsize = pixel_size.to(u.deg).value
        imaging_npixel = 4096 # FIXME: this should be a parameter
        imaging_cellsize = fov / imaging_npixel


        observation_time = sources[0].get_best_observation_time(telescope)
        number_of_timesteps = int(seconds / 7.997) # 7.997 seconds per timestep
        

        observation = Observation(
            start_frequency_hz=start_freq.to(u.Hz).value,
            start_date_and_time=observation_time,
            frequency_increment_hz=delta_freq.to(u.Hz).value,
            length = timedelta(seconds=seconds),
            number_of_time_steps=number_of_timesteps,
            number_of_channels=n_channels,
            phase_centre_ra_deg = x0,
            phase_centre_dec_deg = y0,
        )

        # run a single simulation with the provided configuration 
        simulation = InterferometerSimulation(
            channel_bandwidth_hz=delta_freq.to(u.Hz).value,
            station_type="Gaussian beam",
            gauss_beam_fwhm_deg=fov.to(u.deg).value,
            gauss_ref_freq_hz=frequency.to(u.Hz).value,
            use_gpus=False,
        )


        printlog (log_file, "Starting simulation with params:")
        printlog (log_file, f"\tSources: {len(sources)}")
        for i, src in enumerate(sources):
            printlog (log_file, f"\tSource {i}: {src.to_json()}")
        printlog (log_file, f"\tObservation time: {observation_time.strftime('%Y-%m-%d %H:%M:%S')}")
        printlog (log_file, f"\tPrefix: {uuid_simulation}")
        printlog (log_file, f"\tTelescope: {telescope.name.upper()}")
        printlog (log_file, f"\tFrequency: {frequency.to(u.MHz).value} MHz")
        printlog (log_file, f"\tNumber of channels: {n_channels}")
        printlog (log_file, f"\tDelta frequency: {delta_freq.to(u.MHz)} MHz")
        printlog (log_file, f"\tExposition time: {seconds} seconds")
        printlog (log_file, f"\tPixels: {imaging_npixel}")


        # Get current path
        current_path = os.path.dirname(os.path.realpath(__file__))
        visibility_path = os.path.join(current_path, 'MSDIRS', 
                                       f"{uuid_simulation}.MS") # path to the visibility file; if you use the WSCLEAN backend, the path must be absolute
        printlog (log_file, f"Running simulation with visibility path: {visibility_path}")
        simulation.run_simulation(telescope=telescope, sky=sky, 
                                  observation=observation, 
                                  visibility_path=visibility_path, 
                                  backend=backend)
        printlog (log_file, f"Simulation finished. Visibility file created at {visibility_path}. Recovering visibilities...")
        visibilities = Visibility(visibility_path)


        if backend_name.lower() == "rascil":
            config = RascilDirtyImagerConfig(
                    imaging_npixel=imaging_npixel,
                    imaging_cellsize=imaging_cellsize.to(u.rad).value,
                )
            imager = RascilDirtyImager(config)
            
        else:
            config = OskarDirtyImagerConfig(
                    imaging_phase_centre=SkyCoord(ra=x0, dec=y0, unit="deg"),
                    combine_across_frequencies=True,
                    imaging_npixel=imaging_npixel,
                    imaging_cellsize=imaging_cellsize.to(u.rad).value,
                )
            imager = OskarDirtyImager(config)

        printlog (log_file, f"Creating dirty images...")
        dirty = imager.create_dirty_image(visibilities)
        printlog (log_file, f"Saving dirty image to dirty.png")
        dirty.plot(title=f"Dirty image {backend_name.upper()} ({telescope.name.upper()})", 
                   filename=dirty_path, wcs_enabled=True, xlabel='RA', ylabel='DEC', norm=PowerNorm(0.3))
        printlog (log_file, f"Saved dirty image to {dirty_path}")
        if tofits:
            printlog (log_file, f"Saving dirty image to dirty.fits")
            dirty_fits_path = os.path.join(folder_name, f"{uuid_simulation}_dirty.fits")
            dirty.write_to_file(dirty_fits_path, overwrite=True)
            printlog (log_file, f"Saved dirty fits to {dirty_fits_path}")


        try:
            sources_json = []
            for i, src in enumerate(sources):
                sources_json.append({"ra": src.ra.value, "dec": src.dec.value, "flux": src.I.value, "mute":False, "name":f"source_{i:03}"})
            json_data = {"sources": sources_json, "phase_center": {"ra": x0, "dec": y0}, 
                        "observation_date": observation_time.strftime("%Y-%m-%d %H:%M:%S"), 
                        "frequency": frequency.to(u.Hz).value, "fov": fov.value, "pixel": pixel_size.value, 
                        "telescope": telescope_name, "backend": backend_name, "simulation": uuid_simulation, 
                        "rms": rms, "cleaned": cleaned}
            with open(cfg_path, "w") as f:
                f.write(json.dumps(json_data, indent=4))
            printlog (log_file, f"Saved configuration file to {cfg_path}")
        except Exception as e:
            printlog (log_file, f"Error saving configuration file: {show_exc(e)}")

        # if (cfg is None):
        #     with open(fout.replace(".png", "_sources.json"), "w") as f:
        #         f.write(json.dumps(json_data, indent=4))
        #     printlog (log_file, f"Saved sources to {fout.replace('.png', '_sources.json')}")

        if cleaned:
            printlog (log_file, f"Cleaning the image...{t0}")
            if backend_name.lower() in ["rascil", "all"]:
                try:
                    deconvolved, restored, residual = RascilImageCleaner(
                        RascilImageCleanerConfig(
                            imaging_npixel=imaging_npixel,
                            imaging_cellsize=imaging_cellsize,
                            ingest_vis_nchan=n_channels,
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
                    restored.plot(title=f"Cleaned image (RASCIL) {backend_name.upper()} ({telescope.name.upper()})", filename=cleaned_path, wcs_enabled=True, xlabel='RA', ylabel='DEC')
                    printlog (log_file, f"Saved cleaned image to cleaned.png")
                except Exception as e:
                    printlog (log_file, f"Error cleaning the image: {e}")
                    printlog (log_file, f"Using WSCLEAN instead {t0}")
                    config = WscleanImageCleanerConfig(imaging_npixel=imaging_npixel, imaging_cellsize=pixel_size.to(u.rad).value)
                    cleaner = WscleanImageCleaner(config)
                    if (tofits):
                        print ("Saving fits....")
                        path_fits = os.path.join(folder_name, f"{uuid_simulation}_wscleaned.fits")
                        cleaned = cleaner.create_cleaned_image(visibilities, output_fits_path=path_fits)
                    else:
                        cleaned = cleaner.create_cleaned_image(visibilities)
                    cleaned.plot(title=f"Cleaned image (WSCLEAN) {backend_name.upper()} ({telescope.name.upper()})", filename=cleaned_path, wcs_enabled=True, xlabel='RA', ylabel='DEC')
                    printlog (log_file, f"Saved cleaned image to cleaned.png")
            if backend_name.lower() in ["oskar", "all", "wsclean"]:
                if (backend_name.lower() == "oskar"):
                    printlog (log_file, f"Cleaning not supported for OSKAR, using WSCLEAN")
                config = WscleanImageCleanerConfig(imaging_npixel=imaging_npixel, imaging_cellsize=imaging_cellsize.to(u.rad).value)
                cleaner = WscleanImageCleaner(config)
                if (tofits):
                    path_fits = os.path.join(folder_name, f"{uuid_simulation}_wscleaned.fits")
                    print (f"Saving fits in {path_fits}")
                    cleaned = cleaner.create_cleaned_image(visibilities, output_fits_path=path_fits)
                else:
                    cleaned = cleaner.create_cleaned_image(visibilities)
                
                gamma = 0.3
                cleaned.plot(title=f"Cleaned image (WSCLEAN) {backend_name.upper()} ({telescope.name.upper()})", filename=cleaned_path, wcs_enabled=True, xlabel='RA', ylabel='DEC', norm=PowerNorm(gamma))
                printlog (log_file, f"Saved cleaned image to cleaned.png")
            image3 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "cleaned.png")
        else:
            printlog (log_file, f"Cleaning not requested")
            image3 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "dirty.png")



        image1 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "dirty.png")
        image2 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "sources.png")
        cfg_path = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "cfg.json")

        print("Simulation finished")
        
        


        return JsonResponse({'error': False, 'error_msg': 'Ok', 'image1': image1, 'image2': image2, 'image3': image3, 'uuid': uuid_simulation, 'cfg_path': cfg_path, 'sources':sources_json}) 
    except Exception as e:
        printlog(log_file, f'{show_exc(e)}')
        return JsonResponse({'error': True, 'error_msg': show_exc(e)})

@csrf_exempt
def load(request):
    try:
        uuid_simulation = request.POST.get('uuid', None)
        if uuid_simulation:
            image1 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "dirty.png")
            image2 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "sources.png")
            image3 = os.path.join(settings.STATIC_URL,'simulations', uuid_simulation, "cleaned.png")
            return JsonResponse({'error': False, 'error_msg': 'Ok', 'image1': image1, 'image2': image2, 'image3': image3, 'uuid': uuid_simulation}) 
        else:
            return JsonResponse({'error': True, 'error_msg': 'UUID not found'})
    except Exception as e:
        return JsonResponse({'error': True, 'error_msg': show_exc(e)})

@csrf_exempt
def upload_config(request):
    try:
        if request.method == 'POST':
            # Get the uploaded file
            if ('cfg' not in request.FILES):
                return JsonResponse({'error': True, 'error_msg': 'No file uploaded'})
            
            uploaded_json_file = request.FILES['cfg']
            # Save the file to a temporary location
            json_data = json.loads(uploaded_json_file.read().decode('utf-8'))
            sources = json_data['sources']
            #send json_data to index as post
            json_data = json.dumps(json_data)
            return JsonResponse({'error': False, 'error_msg': 'Ok', 'json_data': json_data})

        elif request.method == 'GET':
            return JsonResponse({'error': True, 'error_msg': 'GET method not allowed'})
    except Exception as e:
        return JsonResponse({'error': True, 'error_msg': show_exc(e)})
    