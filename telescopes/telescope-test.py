import argparse, sys
import matplotlib.pyplot as plt
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
from karabo.simulator_backend import SimulatorBackend

from datetime import datetime
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)

from astropy.wcs import WCS
from karabo.simulation.visibility import *


# from karabo.imaging.imager_base import DirtyImagerConfig
# from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
# from karabo.imaging.imager_rascil import (
#     RascilDirtyImager,
#     RascilDirtyImagerConfig,
#     RascilImageCleaner,
#     RascilImageCleanerConfig,
# )
# from karabo.imaging.imager_wsclean import (
#     WscleanDirtyImager,
#     WscleanImageCleaner,
#     WscleanImageCleanerConfig,
#     create_image_custom_command,
#     )
# from karabo.simulation.visibility import *
# from astropy.coordinates import SkyCoord  # High-level coordinates


# import numpy as np

# from karabo.simulation.interferometer import InterferometerSimulation
# from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.simulation.observation import Observation
# from karabo.util.plotting_util import get_slices

import time

def printlog(msg, t0=None):
    # t0 is the initial timestamp
    if (t0 is not None):
        if t0 == 0:
            print (f"{msg}")
        else:
            date_for_human = datetime.fromtimestamp(t0).strftime('%Y-%m-%d %H:%M:%S')
            print (f"[{date_for_human} + {time.time()-t0:.2f}sec] {msg}")
    else:
        print (f"[{datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}] {msg}")





if __name__ == "__main__":

    t0 = time.time()   

    parser = argparse.ArgumentParser(description="Generate a sky model with sources. If the user does not provide any parameters, the script will generate a random sky model with 10 sources.")
    # parser.add_argument("--telescope", help="Telescope", type=str, default="LOFAR")
    # parser.add_argument("--backend", help="Backend", type=str, default="rascil")

    args = parser.parse_args()

    for _backend in ["rascil", "oskar"]:
        if _backend == "rascil":
            backend = SimulatorBackend.RASCIL
            from karabo.simulation.telescope import RASCILTelescopes
            telescope_types = get_args(RASCILTelescopes)
                

            # if args.telescope in telescope_types:
            #     telescope = Telescope.constructor(args.telescope, backend=backend)
            # else:
            #     printlog (f"Telescope {args.telescope} not found. Loading MEERKAT+ as default. Values accepted: {telescope_types}", t0)
            #     telescope = Telescope.constructor('MEERKAT+', backend=backend)

        else:
            backend = SimulatorBackend.OSKAR
            from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType, OSKAR_TELESCOPE_TO_VERSIONS
            telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)

        for t in telescope_types:
            try:
                printlog (f"Available telescope: {t} - {_backend.upper()}", t0)
                telescope = Telescope.constructor(t, backend=backend)
                # if (args.telescope in teslescope_types):
                #     telescope = Telescope.constructor(args.telescope, backend=backend)
                # else:
                #     printlog (f"Telescope {args.telescope} not found. Loading EXAMPLE as default.", t0)
                #     telescope = Telescope.constructor("EXAMPLE", backend=backend)
                t0 = 0


                printlog (f"Backend loaded: {backend.name.upper()}", t0)
                printlog (f"Telescope loaded: {telescope.name.upper()}", t0)   
                printlog (f"\tTelescope config:", t0)
                printlog (f"\t\tCentre Longitude: {telescope.centre_longitude}", t0)
                printlog (f"\t\tCentre Latitude: {telescope.centre_latitude}", t0)
                printlog (f"\t\tCentre Altitude: {telescope.centre_altitude}", t0)
                telescope.name = f"{telescope.name.upper()} {_backend.upper()}"
                telescope.plot_telescope(file=f"telescope_{telescope.name.upper()}.png".replace(' ','_'))
            except AssertionError as e:
                try:
                    version = None
                    versions = OSKAR_TELESCOPE_TO_VERSIONS[t]
                    for version in versions:
                        printlog (f"Available telescope version: {version}", t0)
                        telescope = Telescope.constructor(t, backend=backend, version=version)
                        t0 = 0


                        printlog (f"Backend loaded: {backend.name.upper()}", t0)
                        printlog (f"Telescope loaded: {telescope.name.upper()}", t0)   
                        printlog (f"\tTelescope config:", t0)
                        printlog (f"\t\tCentre Longitude: {telescope.centre_longitude}", t0)
                        printlog (f"\t\tCentre Latitude: {telescope.centre_latitude}", t0)
                        printlog (f"\t\tCentre Altitude: {telescope.centre_altitude}", t0)
                        telescope.name = f"{telescope.name.upper()} {_backend.upper()} {version}"
                        telescope.plot_telescope(file=f"telescope_{telescope.name.upper()}.png".replace(' ','_'))
                except Exception as e:
                    printlog (f"Error loading version telescope {_backend.upper()}: {e.__class__.__name__} {version} {e}", t0)



            except Exception as e: 
                printlog (f"Error loading telescope {_backend.upper()}: {e.__class__.__name__} {e}", t0)
            

