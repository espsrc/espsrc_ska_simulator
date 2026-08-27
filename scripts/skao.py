import argparse
import multiprocessing as mp
import os
import shutil
import sys
import time
from datetime import datetime, timedelta

import matplotlib
import numpy as np
from matplotlib.colors import PowerNorm

matplotlib.use("Agg")
import warnings

from astropy.utils.exceptions import AstropyDeprecationWarning

warnings.simplefilter("ignore", category=AstropyDeprecationWarning)
warnings.simplefilter("ignore", category=UserWarning)
import astropy.units as u
import karabo
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.imaging.imager_wsclean import (
    WscleanImageCleaner,
    WscleanImageCleanerConfig,
)
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import (
    OSKARTelescopesWithoutVersionType,
    OSKARTelescopesWithVersionType,
    Telescope,
)
from karabo.simulation.visibility import *
from karabo.simulator_backend import SimulatorBackend

# def printlog(fname, *args):
#     # Print to console and log file
#     print (f'[{datetime.now()}]',*args)
#     with open(fname, 'a') as f:
#         print(f'[{datetime.now()}]', *args, file=f)
#         f.flush()
#         os.fsync(f.fileno())
#         f.close()
# def show_exc(exception):
#     exc_type, exc_obj, tb = sys.exc_info()
#     f = tb.tb_frame
#     lineno = tb.tb_lineno
#     filename = f.f_code.co_filename
#     filename_rel = os.path.relpath(filename, os.path.dirname(__file__))
#     app_folder = os.path.basename(os.path.dirname(__file__))
#     return f'EXCEPTION IN ({filename_rel}, LINE {lineno}): {exception} (APP: {app_folder})'
# class Source:
#     def __init__(self, ra, dec, I, Q=0 * u.Jy, U=0 * u.Jy, V=0 * u.Jy, ref_freq=0 * u.Hz, spec_index=0, rot_meas=0 * u.rad/(u.m**2),
#                  major_axis = 0*u.arcsec, minor_axis = 0*u.arcsec, pa=0*u.arcsec,  true_redshift=0, obs_redshift=0, obj_id = None):
#         # Initialize the source with its parameters, checking units
#         list_of_units = [u.deg, u.deg, u.Jy, u.Jy, u.Jy, u.Jy, u.Hz, u.rad/(u.m**2), u.arcsec, u.arcsec, u.arcsec]
#         list_of_values = [ra, dec, I, Q, U, V, ref_freq, rot_meas, major_axis, minor_axis, pa]
#         for i, unit in enumerate(list_of_units):
#             if not isinstance(list_of_values[i], u.Quantity):
#                 list_of_values[i] = list_of_values[i] * unit
#                 print (f"Adding unit {unit} to value {list_of_values[i]}")
#             if list_of_values[i].unit != unit:
#                 try:
#                     list_of_values[i] = list_of_values[i].to(unit)
#                 except u.UnitConversionError:
#                     raise ValueError(f"Value {list_of_values[i]} does not have the correct unit {unit}")
#         self.ra = list_of_values[0]
#         self.dec = list_of_values[1]
#         self.I = list_of_values[2]
#         self.Q = list_of_values[3]
#         self.U = list_of_values[4]
#         self.V = list_of_values[5]
#         self.ref_freq = list_of_values[6]
#         self.spec_index = spec_index
#         self.rot_meas = list_of_values[7]
#         self.major_axis = list_of_values[8]
#         self.minor_axis = list_of_values[9]
#         self.pa = list_of_values[10]
#         self.true_redshift = true_redshift
#         self.obs_redshift = obs_redshift
#         self.obj_id = obj_id
#     @staticmethod
#     def from_name(name):
#         source = acoord.get_icrs_coordinates(name)
#         if source is None:
#             raise ValueError(f"Source {name} not found")
#         return Source(source.ra, source.dec, 1 * u.Jy)
#     def to_json(self):
#         # Convert the source to a JSON serializable dictionary
#         return {
#             "ra": self.ra.to(u.deg).value,
#             "dec": self.dec.to(u.deg).value,
#             "I": self.I.to(u.Jy).value,
#             "Q": self.Q.to(u.Jy).value,
#             "U": self.U.to(u.Jy).value,
#             "V": self.V.to(u.Jy).value,
#             "ref_freq": self.ref_freq.to(u.Hz).value,
#             "spec_index": self.spec_index,
#             "rot_meas": self.rot_meas.value,
#             "major_axis": self.major_axis.to(u.arcsec).value,
#             "minor_axis": self.minor_axis.to(u.arcsec).value,
#             "pa": self.pa.to(u.arcsec).value,
#             "true_redshift": self.true_redshift,
#             "obs_redshift": self.obs_redshift,
#         }
#     def __str__(self):
#         # Return a string representation of the source (only non-zero values)
#         str2print = f"Source(ra={self.ra}, dec={self.dec}, I={self.I}"
#         json_values = self.to_json()
#         for key, value in json_values.items():
#             if key not in ['ra', 'dec', 'I'] and value != 0:
#                 str2print += f", {key}={value}"
#         str2print += ")"
#         return str2print
#     def to_sky_model(self, reduced_form=False):
#         # Convert the source to a SkyModel object
#         if reduced_form:
#             return (self.ra.value, self.dec.value, self.I.value)
#         else:
#             return( self.ra.value, self.dec.value, self.I.value,
#                 self.Q.value, self.U.value, self.V.value,
#                 self.ref_freq.value, self.spec_index,
#                 self.rot_meas.value, self.major_axis.value,
#                 self.minor_axis.value, self.pa.value,
#                 self.true_redshift, self.obs_redshift)
#     def coords(self, frame='icrs'):
#         # Return the coordinates of the source
#         return SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame=frame)
#     def get_best_observation_time(self, telescope:Telescope, date=None):
#         """
#         Returns the local time at which an object with a given RA/Dec culminates (best observation time).
#         Parameters:
#         - ra_hours: Right Ascension in hours (float)
#         - dec_degrees: Declination in degrees (float)
#         - lat_deg: Observer's latitude in degrees (float)
#         - lon_deg: Observer's longitude in degrees (float, positive to the East)
#         - elevation_m: Altitude above sea level (optional)
#         - date: Date as a string 'YYYY-MM-DD' (optional, defaults to today if not provided)
#         - timezone_offset: Time difference relative to UTC (e.g., -6 for CDMX)
#         Returns:
#         - Best time.
#         """
#         if date is None:
#             date = datetime.now().strftime('%Y-%m-%d')
#         coord = SkyCoord(ra=self.ra, dec=self.dec)
#         location = EarthLocation(lat=telescope.centre_latitude * u.deg, lon=telescope.centre_longitude * u.deg, height=telescope.centre_altitude * u.m)
#         midnight = Time(f"{date} 00:00:00") + 12*u.hour  # Mediodía UTC
#         delta = timedelta(minutes=1)
#         best_time = None
#         max_alt = -90
#         for minutes in range(-360, 360):
#             current_time = midnight + minutes * u.minute
#             altaz = coord.transform_to(AltAz(obstime=current_time, location=location))
#             if altaz.alt.deg > max_alt:
#                 max_alt = altaz.alt.deg
#                 best_time = current_time
#         return best_time
from utils import DIAMETERS, Source, printlog, show_exc

if __name__ == "__main__":
    try:
        # Show the PID of the process
        print(f"PID: {os.getpid()}")
        t0 = time.time()
        # Use argparse for define two parameters; first parameter, a list with sources names; second parameters, show or write png
        argparser = argparse.ArgumentParser(description="SKAO simulation script")
        argparser.add_argument(
            "--sources",
            type=str,
            nargs="+",
            help="List of sources to simulate (e.g. 'J2000 00h00m00s +00d00m00s')",
        )
        argparser.add_argument(
            "--prefix", type=str, help="Prefix for filenames", default=None
        )
        argparser.add_argument(
            "--telescope",
            choices=["SKA1LOW", "SKA1MID"],
            type=str,
            help="Telescope to use for the simulation ('SKA1LOW', 'SKA1MID'), default is 'SKA1LOW'",
            default="SKA1LOW",
        )
        argparser.add_argument(
            "--I", type=float, nargs="+", default=None, help="Intensity in Jy"
        )
        argparser.add_argument(
            "--Q", type=float, nargs="+", default=None, help="Q in Jy"
        )
        argparser.add_argument(
            "--U", type=float, nargs="+", default=None, help="U in Jy"
        )
        argparser.add_argument(
            "--V", type=float, nargs="+", default=None, help="V in Jy"
        )
        argparser.add_argument(
            "--ref_freq",
            type=float,
            nargs="+",
            default=None,
            help="Reference frequency in Hz",
        )
        argparser.add_argument(
            "--overwrite", action="store_true", help="Overwrite existing files"
        )
        argparser.add_argument(
            "--freq", type=float, help="Central Freq in MHz", default=1420
        )
        argparser.add_argument(
            "--n_channels", type=int, help="Number of channels", default=4
        )
        argparser.add_argument(
            "--delta_freq", type=float, help="Delta Freq in MHz", default=0.1
        )
        argparser.add_argument(
            "--seconds", type=int, help="Observation Time in seconds", default=1
        )
        argparser.add_argument(
            "--cleaning", action="store_true", help="Use cleaning algorithm"
        )

        argparser.add_argument(
            "--pixels",
            type=int,
            help="Number of pixels in the image",
            default=512,
        )
        argparser.add_argument(
            "--imaging_niter",
            type=int,
            help="Number of iterations for the imager",
            default=1000,
        )

        argparser.add_argument(
            "--fov", type=float, help="Field of view in degrees", default=0
        )
        args = argparser.parse_args()

        frequency = args.freq * u.MHz

        if args.prefix is not None:
            prefix = args.prefix
        else:
            # Prefix in format YYYYMMDD_HHMMSS
            prefix = datetime.now().strftime("%Y%m%d_%H%M")

        log_file = os.path.join(os.path.dirname(__file__), f"{prefix}.log")
        printlog(log_file, "System info")
        printlog(log_file, f"\tKarabo version: {karabo.__version__}")
        try:
            total_memory_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf(
                "SC_PHYS_PAGES"
            )
            total_memory_gb = total_memory_bytes / (1024.0**3)
            printlog(log_file, f"\tTotal RAM: {total_memory_gb:.2f} GB")
        except Exception:
            printlog(log_file, "\tError getting RAM memory (No LINUX OS?)")

        try:
            total_space_in_GB = shutil.disk_usage(
                os.path.expanduser(os.path.dirname(__file__))
            ).total / (1024.0**3)
            printlog(log_file, f"\tTotal disk space: {total_space_in_GB:.2f} GB")
        except Exception:
            printlog(log_file, "\tError getting disk space)")

        printlog(log_file, f"\tNumber of CPU cores: {mp.cpu_count()}")

        if args.fov == 0:
            wavelength = frequency.to(u.m, equivalencies=u.spectral())
            fov = (1.22 * wavelength / DIAMETERS[args.telescope.upper()]) * u.rad
        else:
            fov = (args.fov * u.deg).to(u.rad)

        printlog(log_file, "Starting simulation with params:")
        printlog(log_file, f"\tSources: {args.sources}")
        printlog(log_file, f"\tPrefix: {prefix}")
        printlog(log_file, f"\tTelescope: {args.telescope}")
        printlog(log_file, f"\t\tDiameter: {DIAMETERS['MEERKAT'].value} meters")
        printlog(log_file, f"\t\tFOV: {fov.to(u.deg).value} degrees")
        printlog(log_file, f"\tFrequency: {args.freq} MHz")
        printlog(log_file, f"\tNumber of channels: {args.n_channels}")
        printlog(log_file, f"\tDelta frequency: {args.delta_freq} MHz")
        printlog(log_file, f"\tObservation time: {args.seconds} seconds")
        printlog(log_file, f"\tPixels: {args.pixels}")
        printlog(log_file, f"\tIntensity: {args.I}")
        printlog(log_file, f"\tQ: {args.Q}")
        printlog(log_file, f"\tU: {args.U}")
        printlog(log_file, f"\tV: {args.V}")
        printlog(log_file, f"\tRef freq: {args.ref_freq}")

        sources_names = args.sources
        if not args.sources:
            printlog(log_file, "No sources provided. Exiting.")
            sys.exit(1)

        # Launch the scripts monitor.py (not blocking) with pid of the current process
        import subprocess

        monitor_script = os.path.join(os.path.dirname(__file__), "monitor.py")
        if os.path.exists(monitor_script):
            printlog(
                log_file,
                f"Launching monitor script {monitor_script} with PID {os.getpid()}",
            )
            monitor_proc = subprocess.Popen(
                [sys.executable, monitor_script, str(os.getpid())],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                close_fds=True,
            )

        sources = []
        intensities = args.I
        if intensities is not None:
            if len(intensities) != len(sources_names):
                for i in range(len(sources_names) - len(intensities)):
                    intensities.append(0)
            intensities = [intensity * u.Jy for intensity in intensities]
        else:
            intensities = [
                np.random.randint(1, 100) * u.Jy for _ in range(len(sources_names))
            ]
        for idx, source_name in enumerate(sources_names):
            try:
                source = Source.from_name(source_name)
                source.I = intensities[idx]
                sources.append(source)
            except Exception as e:
                printlog(log_file, show_exc(e))

        try:
            # source = Source.from_name('HCG16')
            skyModel = SkyModel()
            sources_list = []
            for source in sources:
                printlog(log_file, source.to_sky_model(reduced_form=True))
                sources_list.append(source.to_sky_model(reduced_form=True))

            sources_list = np.array(sources_list)

            skyModel.add_point_sources(sources_list)

            printlog(log_file, "Saving sources")

            sources_path = os.path.join(
                os.path.dirname(__file__), f"{prefix}_sources.png"
            )
            skyModel.explore_sky(
                [source.dec.value, source.ra.value],
                filename=sources_path,
                vmin=0,
                vmax=1.05 * np.max(sources_list[:, 2]),
                cfun=np.abs,
            )

            backend = SimulatorBackend.OSKAR
            telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(
                OSKARTelescopesWithoutVersionType
            )
            if args.telescope in telescope_types:
                telescope = Telescope.constructor(args.telescope, backend=backend)
            else:
                printlog(
                    log_file,
                    f"Telescope {args.telescope} not found. Loading SKA1LOW as default. Values accepted: {telescope_types}",
                )
                telescope = Telescope.constructor("SKA1LOW", backend=backend)

            observation_time = sources[0].get_best_observation_time(telescope=telescope)

            start_freq = (
                args.freq * u.MHz - args.n_channels * args.delta_freq / 2 * u.MHz
            )
            delta_freq = args.delta_freq * u.MHz
            n_channels = args.n_channels
            number_of_timesteps = int(args.seconds / 7.997)
            number_of_channels = n_channels

            observation = Observation(
                start_frequency_hz=start_freq.to(u.Hz).value,
                start_date_and_time=observation_time,
                frequency_increment_hz=delta_freq.to(u.Hz).value,
                length=timedelta(seconds=args.seconds),
                number_of_time_steps=number_of_timesteps,
                number_of_channels=n_channels,
                phase_centre_ra_deg=sources[0].ra.value,
                phase_centre_dec_deg=sources[0].dec.value,
            )

            root_path = os.path.dirname(__file__)
            visibility_path = os.path.join(root_path, f"{prefix}_visibilities.MS")
            if os.path.exists(visibility_path):
                if args.overwrite:
                    printlog(
                        log_file,
                        f"Visibility file {visibility_path} already exists. Overwriting it",
                    )
                    shutil.rmtree(visibility_path)
                else:
                    # Ask por confirmation to overwrite
                    print(
                        log_file,
                        f"Visibility file {visibility_path} already exists. Do you want to overwrite it? (y/n)",
                    )
                    answer = input()
                    if answer.lower() != "y":
                        printlog(log_file, "No overwrite, then exiting")
                        sys.exit(0)
                    else:
                        printlog(
                            log_file, f"Overwriting visibility file {visibility_path}"
                        )
                        shutil.rmtree(visibility_path)
            printlog(log_file, f"Creating visibility file {visibility_path}")
            simulation = InterferometerSimulation(
                channel_bandwidth_hz=delta_freq.to(u.Hz).value,
                station_type="Gaussian beam",
                gauss_beam_fwhm_deg=fov.to(u.deg).value,
                gauss_ref_freq_hz=frequency.to(u.Hz).value,
                use_gpus=False,
            )

            simulation.run_simulation(
                telescope=telescope,
                observation=observation,
                sky=skyModel,
                visibility_path=visibility_path,
                backend=backend,
            )
            printlog(log_file, f"Visibilities saved in {visibility_path}")
            printlog(log_file, "Recovering visibilities")
            visibilities = Visibility(visibility_path)

            imaging_cellsize = fov / int(args.pixels)

            config = OskarDirtyImagerConfig(
                imaging_npixel=args.pixels,
                imaging_cellsize=imaging_cellsize.to(u.rad).value,
                combine_across_frequencies=True,
                imaging_phase_centre=sources[0].coords(),
            )
            imager = OskarDirtyImager(config=config)
            dirty_image = imager.create_dirty_image(visibilities)
            dirty_png_path = os.path.join(root_path, f"{prefix}_dirty.png")
            dirty_image.plot(
                title=f"Dirty image {backend.name} ({telescope.name.upper()})",
                wcs_enabled=True,
                xlabel="RA",
                ylabel="DEC",
            )

            dirty_image.plot(
                title=f"Dirty image {backend.name} ({telescope.name.upper()})",
                filename=dirty_png_path,
                wcs_enabled=True,
                xlabel="RA",
                ylabel="DEC",
            )
            printlog(log_file, f"Dirty image (PNG) saved in {dirty_png_path}")

            dirty_fits_path = os.path.join(root_path, f"{prefix}_dirty.fits")
            dirty_image.write_to_file(dirty_fits_path, overwrite=True)
            printlog(log_file, f"Dirty image (FITS) saved in {dirty_fits_path}")

            if args.cleaning:
                printlog(log_file, "Cleaning not supported for OSKAR, using WSCLEAN")
                config = WscleanImageCleanerConfig(
                    imaging_npixel=args.pixels,
                    imaging_cellsize=imaging_cellsize.to(u.rad).value,
                )
                cleaner = WscleanImageCleaner(config)

                path_fits = os.path.join(root_path, f"{prefix}_cleaned.fits")
                cleaned = cleaner.create_cleaned_image(
                    visibilities, output_fits_path=path_fits
                )
                printlog(log_file, f"Cleaned image (FITS) saved in {path_fits}")
                cleaned_path = os.path.join(root_path, f"{prefix}_cleaned.png")
                gamma = 0.3
                cleaned.plot(
                    title=f"Cleaned image (WSCLEAN) {backend.name.upper()} ({telescope.name.upper()})",
                    filename=cleaned_path,
                    wcs_enabled=True,
                    xlabel="RA",
                    ylabel="DEC",
                    norm=PowerNorm(gamma),
                )
                printlog(log_file, "Saved cleaned image to cleaned.png")

            printlog(log_file, "Time ellapsed: ", time.time() - t0)

        except Exception as e:
            printlog(log_file, show_exc(e))
            sys.exit(1)
    except Exception as e:
        printlog(log_file, show_exc(e))
        sys.exit(1)
