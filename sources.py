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



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate a sky model with sources. If the user does not provide any parameters, the script will generate a random sky model with 10 sources.")
    parser.add_argument("--path_cfg", help="Path to the config file", type=str)
    parser.add_argument("--path_out", help="Path to the output file (if this parameter is not provided, then filenames based on the date will be generated)", type=str, default=".png")
    parser.add_argument("--N_srcs", help="Generate N random sources", type=int, default=10)
    parser.add_argument("--rms", help="RMS noise", type=float, default=0.0)
    parser.add_argument("--tofits", help="Save to fits", action="store_true", default=False)
    parser.add_argument("--asksrc", help="Ask for sources", action="store_true", default=False)
    parser.add_argument("--backend", help="Imaging backend", type=str, default="rascil")
    parser.add_argument("--telescope", help="Telescope", type=str, default="LOFAR")

    args = parser.parse_args()

    points = []
    sky = SkyModel()
    sky_data = []

    fout = args.path_out
    if not fout.endswith(".png"):
        fout = fout + ".png"

    if args.path_out == ".png":
        args.path_out = f"sim{time.time():.0f}.png"
        print (f"Output file not specified, saving to {args.path_out}")
        fout = args.path_out


    if args.path_cfg:  # load sources from file
        path_cfg = args.path_cfg
        json_data = json.loads(open(path_cfg).read())
        print (f"Loading sources from file: [{path_cfg}]")
        
        for src in json_data['sources']:
            if "mute" in src:
                if src["mute"]:
                    continue
            sky_data.append([src['ra'], src['dec'], src['flux']])
            points.append((src['ra'], src['dec']))
        sky_data = np.array(sky_data)
        points = [(x, y) for x, y, _ in sky_data]
        x0, y0 = optimal_manhattan_point(points)
        
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
        print (f"Generating {args.N_srcs} random sources")
        limit = 20 * 0.90
        # Random sources [x0, y0, flux]; x0, y0 is the limited between -20, 20. Flux is limited between [0.1, 10] 
        sky_data = np.random.rand(args.N_srcs, 3) * 2*limit - limit
        sky_data[:, 2] = ((sky_data[:, 2] + limit) / (2*limit)) * 9.9 + 0.1
        x0, y0 = 0, 0        
    sky.add_point_sources(sky_data)


    if args.backend == "rascil":
        backend = SimulatorBackend.RASCIL
        from karabo.simulation.telescope import RASCILTelescopes
        telescope_types = get_args(RASCILTelescopes)
        if args.telescope in telescope_types:
            telescope = Telescope.constructor(args.telescope)
        else:
            print (f"Telescope {args.telescope} not found. Loading MEERKAT+ as default.")
            telescope = Telescope.constructor("MEERKAT+")

    else:
        backend = SimulatorBackend.OSKAR
        from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType
        teslescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
        if (args.telescope in teslescope_types):
            telescope = Telescope.constructor(args.telescope)
        else:
            print (f"Telescope {args.telescope} not found. Loading EXAMPLE as default.")
            telescope = Telescope.constructor("EXAMPLE")

    print (f"Telescope loaded: {telescope.name}")

   
    # overwrite or set any of the implemented configuration values
    telescope.centre_longitude = 0.0
    telescope.centre_latitude = 0.0
    telescope.centre_altitude = 0.0
    simulation = InterferometerSimulation()

    # create new observational settings 
    observation = Observation(
        start_frequency_hz=1e6,
        start_date_and_time=datetime(2024, 3, 15, 10, 46, 0),
        phase_centre_ra_deg = x0,
        phase_centre_dec_deg = y0,
    )

    # run a single simulation with the provided configuration
    visibility_path = "./aux.MS" # path to the visibility file
    simulation.run_simulation(telescope, sky, observation, visibility_path=visibility_path)

    visibilities = Visibility(visibility_path)

    imaging_npixel = 2048 
    imaging_cellsize = 3.878509448876288e-03 * 0.1

    config = RascilDirtyImagerConfig(
            imaging_npixel=imaging_npixel,
            imaging_cellsize=imaging_cellsize,
        )
    dirty = RascilDirtyImager(config).create_dirty_image(visibilities)

    if args.rms > 0:
        print (f"Adding noise to the image with RMS: {args.rms}")
        noise_matrix = np.random.normal(0, args.rms, dirty.data[0].shape)
        dirty.data[0] += noise_matrix


    dirty.plot(title="Dirty image RASCIL", filename=fout, wcs_enabled=True)
    print (f"Sky Model with {len(sky_data)} sources")
    for src in sky_data:
        print (f"Source: {src[0]}, {src[1]}, {src[2]}")
    print (f"Optimal phase center: {x0}, {y0}")

    wcs = WCS(dirty.header)
    slices = get_slices(wcs=wcs)
    fig, ax = plt.subplots(
                 subplot_kw=dict(projection=wcs, slices=slices)
            )
    
    sky_data = np.array(sky_data)

    model =  np.zeros((imaging_npixel, imaging_npixel))
    for src in sky_data:
        x, y, flux = src
        [x, y, _, _] = wcs.all_world2pix(x, y, 0, 0, 0)
        model[int(x), int(y)] = flux
    im = ax.imshow(model, origin='lower', cmap='jet')
    
    ax.set_title("Sources position")
    max_flux = np.max(sky_data[:, 2])
    cmap = plt.get_cmap('jet')

    for src in sky_data:
        x, y, flux = src
        coord = SkyCoord(ra=x, dec=y, unit="deg")
        flux_color = cmap(flux/max_flux)
        [x, y, _, _] = wcs.all_world2pix(x, y, 0, 0, 0)
        circulo = plt.Circle((int(x), int(y)), flux*10, color=flux_color, fill=True, linewidth=2)
        ax.add_patch(circulo)
    fig.colorbar(im)
    fig.savefig(fout.replace(".png", "_sources.png"))

    # Convert sky_data to json
    sky_data = sky_data.tolist()
    sources = []
    for i, src in enumerate(sky_data):
        sources.append({"ra": src[0], "dec": src[1], "flux": src[2], "mute":False, "name":f"source_{i:03}"})
    json_data = {"sources": sources}
    with open(fout.replace(".png", "_sources.json"), "w") as f:
        f.write(json.dumps(json_data, indent=4))
    print (f"Saved sources to {fout.replace('.png', '_sources.json')}")

    if args.tofits:
        dirty.write_to_file(fout.replace(".png", ".fits"))
        print (f"Saved image to {fout.replace('.png', '.fits')}")

