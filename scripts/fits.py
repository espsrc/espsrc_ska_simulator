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
from astropy.wcs import WCS
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
            choices=["SKA1LOW", "SKA1MID", "MeerKAT"],
            type=str,
            help="Telescope to use for the simulation ('SKA1LOW', 'SKA1MID'), default is 'SKA1LOW'",
            default="SKA1LOW",
        )
        argparser.add_argument(
            "--I", type=float, default=None, help="Total Intensity in Jy"
        )
        argparser.add_argument("--Q", type=float, default=None, help="Q")
        argparser.add_argument("--U", type=float, default=None, help="U")
        argparser.add_argument("--V", type=float, default=None, help="V")
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
            "--freq", type=float, help="Central Freq in MHz", default=200
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
        argparser.add_argument(
            "--rms", action="store_true", default=False, help="Enable RMS calculation"
        )
        # Add noise_level argument
        argparser.add_argument(
            "--rms_start",
            type=float,
            default=0,
            help="Noise level start in Jy for RMS calculation",
        )
        # Add noise_level_var argument
        argparser.add_argument(
            "--rms_end",
            type=float,
            default=0,
            help="Noise level end in Jy for RMS calculation",
        )
        args = argparser.parse_args()

        frequency = args.freq * u.MHz

        if args.prefix is not None:
            prefix = args.prefix
        else:
            # Prefix in format YYYYMMDD_HHMMSS
            prefix = datetime.now().strftime("%Y%m%d_%H%M")

        log_file = os.path.join(os.path.dirname(__file__), f"{prefix}.log")
        printlog(log_file, f"Command line: {' '.join(sys.argv)}")
        printlog(log_file, "System info")
        # Print command line
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
        if args.fov == 0:
            # Convert frequency to wavelength
            wavelength = frequency.to(u.m, equivalencies=u.spectral())
            fov = (1.22 * wavelength / DIAMETERS[args.telescope.upper()]) * u.rad
        else:
            fov = (args.fov * u.deg).to(u.rad)

        printlog(log_file, "Starting simulation with params:")
        printlog(log_file, f"\tSources: {args.sources}")
        printlog(log_file, f"\tPrefix: {prefix}")
        printlog(log_file, f"\tTelescope: {args.telescope}")
        printlog(
            log_file, f"\t\tDiameter: {DIAMETERS[args.telescope.upper()].value} meters"
        )
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
        printlog(log_file, f"\tCleaning: {args.cleaning}")
        printlog(log_file, f"\tRMS: {args.rms}")

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
                [
                    sys.executable,
                    monitor_script,
                    str(os.getpid()),
                    f"--csv={prefix}_monitor.log",
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                close_fds=True,
            )

        try:
            source_ref = Source.from_name("HCG16")
            sources_path = os.path.join(
                os.path.dirname(__file__), f"{prefix}_sources.png"
            )
            survey_file = os.path.join(
                os.path.dirname(__file__), "skymodel_incl0_stokes.fits"
            )
            from astropy.io import fits

            if os.path.exists(survey_file):
                fits_data = fits.open(survey_file)
                fits_header = fits_data[0].header
                fits_data = fits_data[0].data
                img_pixels = int(fits_data.shape[2])
                fits_wcs = WCS(fits_header)
                sky_wcs = WCS(naxis=4)  # RA, DEC, Intensities, STOKES
                sky_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ", "STOKES"]
                sky_wcs.wcs.crpix = [
                    fits_data.shape[2] // 2,
                    fits_data.shape[3] // 2,
                    0,
                    0,
                ]
                sky_wcs.wcs.crval = [
                    source_ref.ra.to(u.deg).value,
                    source_ref.dec.to(u.deg).value,
                    frequency.to(u.Hz).value,
                    1.0,
                ]
                sky_wcs.wcs.cdelt = [
                    fov.to(u.deg).value / img_pixels,
                    fov.to(u.deg).value / img_pixels,
                    1.0,
                    1.0,
                ]
                sky_wcs.wcs.cunit = ["deg", "deg", "Hz", ""]

                fluxes = fits_data[0, 0, :, :]
                if args.I is not None and args.I > 0:
                    fluxes = (
                        fluxes / np.max(fluxes) * args.I
                    )  # Normalize to max intensity

                sources = []
                total_pixels = fluxes.size
                skyModel = SkyModel(wcs=sky_wcs)

                fluxes_nonzero = np.nonzero(fluxes)
                indices = np.array(fluxes_nonzero).T

                progress = 0
                total_pixels = indices.shape[0]
                progress_to_print = np.linspace(0, total_pixels, 11, dtype=int)
                max_flux = np.max(fluxes)
                t0 = time.time()
                printlog(log_file, "Starting conversion...")
                ra_list = []
                dec_list = []
                flux_list = []
                sum_weights = np.sum(fluxes)
                if args.I is not None and args.I > 0:
                    fluxes = fluxes / sum_weights * args.I
                elif args.I is None or args.I <= 0:
                    fluxes = fluxes / sum_weights * 1  # Normalize to 1 Jy

                for x, y in indices:
                    world = sky_wcs.pixel_to_world(x, y, 0, 0)
                    skycoord, freq, _ = world
                    intensity = fluxes[x, y] * u.Jy
                    ra_list.append(skycoord.ra.value)
                    dec_list.append(skycoord.dec.value)
                    flux_list.append(intensity.value)
                    progress += 1
                    print(
                        f"Progress: {progress:5.0f}/{total_pixels} ({progress / total_pixels * 100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds",
                        end="\r",
                    )

                    if progress in progress_to_print:
                        printlog(
                            log_file,
                            f"Progress: {progress:5.0f}/{total_pixels} ({progress / total_pixels * 100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds",
                        )

                # np_samples = np.vstack((np.array(ra_list), np.array(dec_list))).transpose()
                np_samples = np.vstack((np.array(ra_list), np.array(dec_list)))
                np_fluxes = np.reshape(np.array(flux_list), (len(flux_list), 1))
                sky_array = np.hstack((np_samples, np_fluxes))
                skyModel = SkyModel(sky_array, wcs=sky_wcs)

                # pickle_path = os.path.join(os.path.dirname(__file__), f'{prefix}_sky_model.pkl')
                # with open(pickle_path, 'wb') as f:
                #     pickle.dump(skyModel.sources, f)
                # printlog (log_file, f"Sky model saved in {pickle_path}")

            printlog(log_file, "Saving sources")
            skyModel.explore_sky(
                [source_ref.dec.value, source_ref.ra.value],
                filename=sources_path,
                cfun=np.abs,
            )
            printlog(log_file, f"Sources saved in {sources_path}. Finished.")
            observation_time = source_ref.get_best_observation_time(telescope=telescope)
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
                phase_centre_ra_deg=source_ref.ra.value,
                phase_centre_dec_deg=source_ref.dec.value,
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

            simulation_params = {
                "channel_bandwidth_hz": delta_freq.to(u.Hz).value,
                #'station_type': "Gaussian beam",
                "sation_type": "Aperture array",
                "gauss_beam_fwhm_deg": fov.to(u.deg).value,
                "gauss_ref_freq_hz": frequency.to(u.Hz).value,
                "use_gpus": False,
            }

            if args.rms:
                simulation_params["noise_enable"] = True
                simulation_params["noise_freq"] = "Observation settings"
                simulation_params["noise_rms_start"] = args.rms_start
                simulation_params["noise_rms_end"] = args.rms_end

            simulation = InterferometerSimulation(**simulation_params)

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
                imaging_phase_centre=source_ref.coords(),
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
        print(show_exc(e))
        sys.exit(1)
