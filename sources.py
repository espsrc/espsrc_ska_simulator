import argparse
import json, sys
import matplotlib.pyplot as plt
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

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

import time

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

from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz
import astropy.units as u

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

# Example usage
longitude = 116.76444824  # degrees East
latitude = -26.82472208   # degrees South
observation_time_local = "2025-03-06 23:59:00"  # Local time at the telescope

result = calculate_lst(longitude, latitude, observation_time_local)
print(f"Local Sidereal Time (LST): {result['LST']} h")
print(f"RA of the zenith: {result['RA_zenith']:.4f}°")
print(f"Dec of the zenith: {result['Dec_zenith']:.4f}°")



if __name__ == "__main__":

    t0 = time.time()   

    parser = argparse.ArgumentParser(description="Generate a sky model with sources. If the user does not provide any parameters, the script will generate a random sky model with 10 sources.")
    parser.add_argument("--cfg", help="Path to the config file", type=str, default=None)
    parser.add_argument("--out", help="Path to the output file (if this parameter is not provided, then filenames based on the date will be generated)", type=str, default=".png")
    parser.add_argument("--N_srcs", help="Generate N random sources", type=int, default=10)
    parser.add_argument("--rms", help="RMS noise", type=float, default=0.0)
    parser.add_argument("--tofits", help="Save to fits", action="store_true", default=False)
    parser.add_argument("--asksrc", help="Ask for sources", action="store_true", default=False)
    parser.add_argument("--backend", help="Imaging backend", type=str, default="rascil")
    parser.add_argument("--telescope", help="Telescope", type=str, default="LOFAR")
    parser.add_argument("--maxflux", help="Max flux", type=float, default=10.)
    parser.add_argument("--freq", help="Frequency in MHz (default = 1 MHz)", type=float, default=1e6) 
    parser.add_argument("--fov", help="Field of view in degrees (default = 1 degree)", type=float, default=1.0) 
    parser.add_argument("--pixel", help="Size of pixel in arcsec (default = 1)", type=float, default=1.0) 
    parser.add_argument("--interactive", help="Interactive mode", action="store_true", default=False)


    args = parser.parse_args()






    fout = args.out
    if not fout.endswith(".png"):
        fout = fout + ".png"

    if args.out == ".png":
        args.out = f"sim{time.time():.0f}.png"
        printlog (f"Output file not specified, saving to {args.out}", t0)
        fout = args.out

    if args.backend == "rascil":
        backend = SimulatorBackend.RASCIL
        from karabo.simulation.telescope import RASCILTelescopes
        telescope_types = get_args(RASCILTelescopes)
        if args.telescope in telescope_types:
            telescope = Telescope.constructor(args.telescope, backend=backend)
        else:
            printlog (f"Telescope {args.telescope} not found. Loading LOFAR as default. Values accepted: {telescope_types}", t0)
            telescope = Telescope.constructor('LOFAR', backend=backend)

    else:
        backend = SimulatorBackend.OSKAR
        from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType
        telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
        if (args.telescope in telescope_types):
            telescope = Telescope.constructor(args.telescope, backend=backend)
        else:
            printlog (f"Telescope {args.telescope} not found. Loading LOFAR as default. Values accepted: {telescope_types}", t0)
            telescope = Telescope.constructor("LOFAR", backend=backend)
    printlog (f"Telescope loaded: {telescope.name.upper()}", t0)

    observation_time_local = datetime.now()

    points = []
    sky_data = []
    printlog(f"Starting the simulation. Command line arguments: {' '.join(sys.argv)}", t0)

    if args.cfg:  # load sources from file
        path_cfg = args.cfg
        json_data = json.loads(open(path_cfg).read())
        printlog (f"Loading sources from file: [{path_cfg}]", t0)
        
        for src in json_data['sources']:
            if "mute" in src:
                if src["mute"]:
                    continue
            sky_data.append([src['ra'], src['dec'], src['flux']])
            points.append((src['ra'], src['dec']))
        sky_data = np.array(sky_data)
        points = [(x, y) for x, y, _ in sky_data]
        x0, y0 = json_data['phase_center']['ra'], json_data['phase_center']['dec']
        
    elif (args.asksrc):
        print ("Enter the sources in the format: ra dec flux")
        print ("Enter 'q' to finish")
        while True:
            src = input()
            if src == "q":
                break
            try:
                ra, dec, flux = src.split()
                ra, dec, flux = float(ra), float(dec), float(flux)
                sky_data.append([ra, dec, flux])
                points.append((ra, dec))
            except:
                print ("Invalid format")
                continue
        sky_data = np.array(sky_data)
        points = [(x, y) for x, y, _ in sky_data]
        x0, y0 = optimal_manhattan_point(points)
    else:        # create a simple sky model with three point sources
        printlog (f"Generating {args.N_srcs} random sources", t0)
        limit = args.fov/2
        result = calculate_lst(telescope.centre_longitude, telescope.centre_latitude, observation_time_local.strftime("%Y-%m-%d %H:%M:%S"))
        x0 = result["RA_zenith"]
        y0 = result["Dec_zenith"]
        # Random sources [x0, y0, flux]; x0, y0 is the limited between -20, 20. Flux is limited between [0.1, 10] 
        sky_data = np.random.rand(args.N_srcs, 3) * 2*limit - limit
        sky_data[:, 2] = ((sky_data[:, 2] + limit) / (2*limit)) * (args.maxflux - 0.1) + 0.1
        sky_data[:, 0] += x0
        sky_data[:, 1] += y0

    observation = Observation(
        start_frequency_hz=1e6,
        start_date_and_time=observation_time_local,
        phase_centre_ra_deg = x0,
        phase_centre_dec_deg = y0,
    )

    imaging_cellsize = (args.pixel *  u.arcsec).to(u.deg).value
    imaging_npixel = 2048

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
    sky.explore_sky([x0, y0],wcs=wcs, filename=fout.replace(".png", "_sources.png"), xlabel='RA', ylabel='DEC')

    # run a single simulation with the provided configuration 
    simulation = InterferometerSimulation()
    visibility_path = "./aux.MS" # path to the visibility file
    simulation.run_simulation(telescope, sky, observation, visibility_path=visibility_path)
    visibilities = Visibility(visibility_path)


    # Create the dirty image
    imaging_cellsize = (args.pixel *  u.arcsec).to(u.deg).value
    imaging_npixel = 2048

    if args.backend == "rascil":
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
    if args.rms > 0:
        printlog (f"Adding noise to the image with RMS: {args.rms}", t0)
        noise_matrix = np.random.normal(0, args.rms, dirty.data[0].shape)
        dirty.data[0] += noise_matrix


    dirty.plot(title=f"Dirty image {args.backend.upper()} ({telescope.name.upper()})", filename=fout, wcs_enabled=True, xlabel='RA', ylabel='DEC')

    printlog (f"Sky Model with {len(sky_data)} sources", t0)
    printlog (f"Optimal phase center: {x0}, {y0}", t0)

    sky_data = np.array(sky_data)

    # Convert sky_data to json
    sky_data = sky_data.tolist()
    sources = []
    for i, src in enumerate(sky_data):
        sources.append({"ra": src[0], "dec": src[1], "flux": src[2], "mute":False, "name":f"source_{i:03}"})
    json_data = {"sources": sources, "phase_center": {"ra": x0, "dec": y0}, "observation_date": observation_time_local.strftime("%Y-%m-%d %H:%M:%S")}
    if (args.cfg is None):
        with open(fout.replace(".png", "_sources.json"), "w") as f:
            f.write(json.dumps(json_data, indent=4))
        printlog (f"Saved sources to {fout.replace('.png', '_sources.json')}", t0)

    if args.tofits:
        dirty.write_to_file(fout.replace(".png", ".fits"), overwrite=True)
        printlog (f"Saved image to {fout.replace('.png', '.fits')}", t0)
        
    printlog ("Cleaning the image...", t0)
    if args.backend == "rascil":
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
        ).create_cleaned_image_variants(visibilities)
        restored.plot(title=f"Cleaned image {args.backend.upper()} ({telescope.name.upper()})", filename=fout.replace(".png", "_cleaned.png"), wcs_enabled=True, xlabel='RA', ylabel='DEC')
        printlog (f"Saved cleaned image to {fout.replace('.png', '_cleaned.png')}", t0)
    else:
        printlog ("Cleaning not supported for OSKAR, using WSCLEAN", t0)
        config = WscleanImageCleanerConfig(
            imaging_npixel=imaging_npixel,
            imaging_cellsize=imaging_cellsize,
            clean_threshold=0.12e-3,
            clean_niter=1000,
            clean_gain=0.1,
            clean_weighting="briggs",
            clean_robust=0.5,
            clean_psf_support=640,
            clean_restored_output="integrated",
        )
        cleaner = WscleanImageCleaner(config)
        cleaned = cleaner.create_cleaned_image_variants(visibilities)
        cleaned.plot(title=f"Cleaned image {args.backend.upper()} ({telescope.name.upper()})", filename=fout.replace(".png", "_cleaned.png"), wcs_enabled=True, xlabel='RA', ylabel='DEC')
        printlog (f"Saved cleaned image to {fout.replace('.png', '_cleaned.png')}", t0)

