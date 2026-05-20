"""pipeline.py — end-to-end interferometric simulation orchestrator."""

from __future__ import annotations

import glob
import json
import os
import pickle
import shutil
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import UnitBase
from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
from karabo.simulation.sky_model import SkyPrefixMapping, SkySourcesUnits
from karabo.simulation.telescope import Telescope
from karabo.simulator_backend import SimulatorBackend
from loguru import logger

from .config import SimConfig
from .imaging import run_dirty_imaging, run_wsclean_imaging
from .manifest import RunContext, create_run_context
from .sky import SkyModel, Source
from .utils import get_diameter

# --------------------------------------------------------------------------- #
# workdir + logging
# --------------------------------------------------------------------------- #


# (setup_workdir moved to create_run_context in manifest.py)


# --------------------------------------------------------------------------- #
# telescope
# --------------------------------------------------------------------------- #


def build_telescope(ctx: RunContext):
    """return a Karabo Telescope instance."""
    config = ctx.config
    kwargs: dict = {"backend": SimulatorBackend.OSKAR}
    if config.telescope_version is not None:
        kwargs["version"] = config.telescope_version
        logger.info(f"Telescope {config.telescope}  version={config.telescope_version}")
    else:
        logger.info(f"Telescope {config.telescope}  (no version)")
    telescope = Telescope.constructor(config.telescope, **kwargs)
    ctx.add_milestone(
        "telescope_built",
        "completed",
        details={"name": config.telescope, "version": config.telescope_version},
    )
    return telescope


# --------------------------------------------------------------------------- #
# field of view
# --------------------------------------------------------------------------- #


def compute_fov(config: SimConfig, frequency: u.Quantity) -> u.Quantity:
    """return FoV in radians.  If fov_deg is set, use it; else diffraction limit."""
    if config.imaging.fov_deg is not None:
        return (config.imaging.fov_deg * u.deg).to(u.rad)
    wavelength = frequency.to(u.m, equivalencies=u.spectral())
    diameter = get_diameter(config.telescope.upper())
    fov = (1.25 * wavelength / diameter) * u.rad
    logger.debug(f"Computed FOV is: {fov} radians")
    return fov


# --------------------------------------------------------------------------- #
# phase centre
# --------------------------------------------------------------------------- #


def parse_center(center_str: Optional[str], fallback: SkyCoord) -> SkyCoord:
    if center_str is None:
        return fallback
    try:
        coords_str = center_str.replace(",", " ").replace(":", " ")
        return SkyCoord(coords_str, unit=(u.hourangle, u.deg))
    except Exception:
        return fallback


# --------------------------------------------------------------------------- #
# sky model — file loaders
# --------------------------------------------------------------------------- #


def _load_sky_from_file(
    fpath: str,
    column_mapping: str = "0,1,2,3,4,5,6,7,8,9,10,11,12",
    scale_I: float = 1.0,
    ref_freq_hz: Optional[float] = None,
    frequency: Optional[u.Quantity] = None,
) -> SkyModel:
    """Load SkyModel from pickle, fits, json, or karabo.mod."""
    ext = os.path.splitext(fpath)[-1].lower()

    # pickle / karabo model
    if ext in (".pkl", ".pickle", ".kmod", ".karabo.mod"):
        with open(fpath, "rb") as fh:
            sky_model = pickle.load(fh)
        logger.info(f"Loaded pickle model from {fpath}")
        return sky_model

    # json catalogue
    if ext == ".json":
        with open(fpath, "r") as fh:
            data = json.load(fh)
        sources: List[Source] = []
        for item in data:
            src = Source.from_json(item)
            if scale_I != 1.0:
                src.I *= scale_I
            if src.ref_freq == 0:
                hz = (
                    ref_freq_hz
                    if ref_freq_hz is not None
                    else (frequency.to(u.Hz).value if frequency is not None else 0)
                )
                src.ref_freq = hz * u.Hz
            sources.append(src)
        if not sources:
            raise ValueError(f"No sources found in JSON {fpath}")
        sky_model = SkyModel()
        arr = np.array([s.to_sky_model(reduced_form=True) for s in sources])
        sky_model.add_point_sources(arr)
        sky_model.get_center()
        logger.info(f"Loaded {len(sources)} sources from JSON {fpath}")
        return sky_model

    # fits table or image
    if ext in (".fits", ".fit"):
        return _load_sky_from_fits(
            fpath, column_mapping, scale_I, ref_freq_hz, frequency
        )

    raise ValueError(f"Unsupported sky-file extension: {ext}")


def _load_sky_from_fits(
    fpath: str,
    column_mapping: str,
    scale_I: float,
    ref_freq_hz: Optional[float],
    frequency: Optional[u.Quantity],
) -> SkyModel:
    """try Karabo's get_sky_model_from_fits; fallback to our Source.from_table_in_fits."""
    cols_mapping = [int(i) for i in column_mapping.split(",")]

    with fits.open(fpath) as hdul:
        hdu1 = hdul[1]
        unit_mapping: Dict[str, UnitBase] = {}
        for col in hdu1.columns:
            mapped = mapping_unit(col.unit)
            unit_mapping[col.unit] = (
                u.Unit(mapped) if mapped else u.dimensionless_unscaled
            )

        prefix_mapping = SkyPrefixMapping(
            ra=hdu1.columns.names[cols_mapping[1]],
            dec=hdu1.columns.names[cols_mapping[2]],
            stokes_i=hdu1.columns.names[cols_mapping[3]],
            stokes_q=hdu1.columns.names[cols_mapping[4]]
            if cols_mapping[4] > -1
            else None,
            stokes_u=hdu1.columns.names[cols_mapping[5]]
            if cols_mapping[5] > -1
            else None,
            stokes_v=hdu1.columns.names[cols_mapping[6]]
            if cols_mapping[6] > -1
            else None,
            spectral_index=hdu1.columns.names[cols_mapping[7]]
            if cols_mapping[7] > -1
            else None,
            ref_freq=hdu1.columns.names[cols_mapping[8]]
            if cols_mapping[8] > -1
            else None,
            rm=hdu1.columns.names[cols_mapping[9]] if cols_mapping[9] > -1 else None,
            major=hdu1.columns.names[cols_mapping[10]]
            if cols_mapping[10] > -1
            else None,
            minor=hdu1.columns.names[cols_mapping[11]]
            if cols_mapping[11] > -1
            else None,
            pa=hdu1.columns.names[cols_mapping[12]] if cols_mapping[12] > -1 else None,
            id=hdu1.columns.names[cols_mapping[0]] if cols_mapping[0] > -1 else None,
        )

        units_sources = SkySourcesUnits(
            stokes_i=u.Jy / u.beam,
            stokes_q=u.Jy / u.beam,
            stokes_u=u.Jy / u.beam,
            stokes_v=u.Jy / u.beam,
            ref_freq=u.MHz,
            major=u.arcsec,
            minor=u.arcsec,
            pa=u.deg,
            rm=u.rad / u.m**2,
        )

    # first attempt with beam units
    try:
        sky_model = SkyModel.get_sky_model_from_fits(
            fits_file=fpath,
            prefix_mapping=prefix_mapping,
            unit_mapping=unit_mapping,
            units_sources=units_sources,
            min_freq=None,
            max_freq=None,
            encoded_freq=None,
            memmap=False,
        )
        logger.info(f"Loaded FITS via Karabo: {fpath}")
        return sky_model
    except u.core.UnitConversionError as exc:
        logger.error(f"Beam-unit conversion failed ({exc}); retrying without beam.")

    # retry without /beam
    units_sources = SkySourcesUnits(
        stokes_i=u.Jy,
        stokes_q=u.Jy,
        stokes_u=u.Jy,
        stokes_v=u.Jy,
        ref_freq=u.MHz,
        major=u.arcsec,
        minor=u.arcsec,
        pa=u.deg,
        rm=u.rad / u.m**2,
    )
    sky_model = SkyModel.get_sky_model_from_fits(
        fits_file=fpath,
        prefix_mapping=prefix_mapping,
        unit_mapping=unit_mapping,
        units_sources=units_sources,
        min_freq=None,
        max_freq=None,
        encoded_freq=None,
        memmap=False,
    )
    logger.info(f"Loaded FITS via Karabo (no-beam fallback): {fpath}")
    return sky_model


# --------------------------------------------------------------------------- #
# sky model — high-level builder
# --------------------------------------------------------------------------- #


def build_sky_model(
    ctx: RunContext,
    fov: u.Quantity,
) -> tuple[SkyModel, SkyCoord]:
    """Return (sky_model, center)."""
    config = ctx.config

    # 1) file path given?
    if config.sky_file is not None:
        sky_model = _load_sky_from_file(
            str(ctx.sky_file_resolved),
            column_mapping=config.column_mapping or "0,1,2,3,4,5,6,7,8,9,10,11,12",
            scale_I=config.scale_I,
            ref_freq_hz=(config.ref_freq_hz[0] if config.ref_freq_hz else None),
            frequency=config.observation.freq_mhz * u.MHz,
        )
        center = sky_model.get_center()
        n_srcs = len(sky_model.sources) if hasattr(sky_model, "sources") else None
        ctx.add_milestone(
            "sky_model_loaded",
            "completed",
            details={
                "path": str(ctx.sky_file_resolved),
                "format": config.sky_format,
                "n_sources": n_srcs,
            },
        )
        return sky_model, center

    # 2) built-in catalogue
    if config.catalogue > 0:
        if config.catalogue == 1:
            logger.info("Loading MIGHTEE catalogue")
            sky_model = SkyModel.get_MIGHTEE_Sky()
            fmt = "MIGHTEE"
        elif config.catalogue == 2:
            logger.info("Loading GLEAM catalogue")
            sky_model = SkyModel.get_GLEAM_Sky()
            fmt = "GLEAM"
        elif config.catalogue == 3:
            skamid_path = Path("SKAMid_B1_8h_v3.fits").resolve()
            if skamid_path.exists():
                logger.info(f"Loading SKAMid catalogue {skamid_path}")
                sky_model = SkyModel.get_sky_model_from_fits(fits_file=str(skamid_path))
                fmt = "SKAMid"
            else:
                logger.info(f"SKAMid catalogue not found at {skamid_path}")
                raise FileNotFoundError(str(skamid_path))
        else:
            raise ValueError(
                f"Catalogue {config.catalogue} not available (1=MIGHTEE, 2=GLEAM, 3=SKAMid)"
            )
        center = sky_model.get_center()
        n_srcs = len(sky_model.sources) if hasattr(sky_model, "sources") else None
        ctx.add_milestone(
            "sky_model_loaded",
            "completed",
            details={"format": fmt, "n_sources": n_srcs},
        )
        return sky_model, center

    # 3) random sources around a reference position
    logger.info("Generating random sources")
    source_ref = Source.from_name("HCG16")
    intensities = [i * u.Jy for i in config.I]
    n_sources = len(intensities)
    sources: List[Source] = []
    for idx in range(n_sources):
        if idx == 0:
            src = source_ref
            src.I = intensities[idx]
        else:
            x_coord = np.random.uniform(-fov.value / 2, fov.value / 2) * 0.8 * u.rad
            y_coord = np.random.uniform(-fov.value / 2, fov.value / 2) * 0.8 * u.rad
            src = Source(
                source_ref.ra + x_coord, source_ref.dec + y_coord, intensities[idx]
            )
        sources.append(src)

    sky_model = SkyModel()
    arr = np.array([s.to_sky_model(reduced_form=True) for s in sources])
    sky_model.add_point_sources(arr)
    center = sky_model.get_center()
    ctx.add_milestone(
        "sky_model_loaded",
        "completed",
        details={"format": "random", "n_sources": n_sources, "reference": "HCG16"},
    )
    return sky_model, center


# --------------------------------------------------------------------------- #
# observation
# --------------------------------------------------------------------------- #


def build_observation(
    ctx: RunContext,
    center: SkyCoord,
    telescope,
) -> tuple:
    """Return (observation, frequency, bandwidth, n_channels, delta_freq, start_freq)."""
    config = ctx.config
    obs = config.observation
    freq = obs.freq_mhz * u.MHz
    bw_mhz = obs.bandwidth_mhz
    n_channels = obs.n_channels
    df_mhz = obs.delta_freq_mhz
    bandwidth = bw_mhz * u.MHz
    delta_freq = df_mhz * u.MHz
    seconds = config.observation.seconds

    start_freq = freq - n_channels * delta_freq / 2

    # best observation time (culmination)
    obs_time = source_ref_get_best_observation_time(center, telescope)
    n_timesteps = max(1, int(seconds / 7.997))

    observation = Observation(
        start_frequency_hz=start_freq.to(u.Hz).value,
        start_date_and_time=obs_time,
        frequency_increment_hz=delta_freq.to(u.Hz).value,
        length=timedelta(seconds=seconds),
        number_of_time_steps=n_timesteps,
        number_of_channels=n_channels,
        phase_centre_ra_deg=center.ra.to(u.deg).value,
        phase_centre_dec_deg=center.dec.to(u.deg).value,
    )

    ctx.add_milestone(
        "observation_configured",
        "completed",
        details={
            "freq_mhz": obs.freq_mhz,
            "bandwidth_mhz": bw_mhz,
            "n_channels": n_channels,
            "seconds": seconds,
        },
    )
    return observation, freq, bandwidth, n_channels, delta_freq, start_freq


def source_ref_get_best_observation_time(center: SkyCoord, telescope):
    """Wrapper around Source.get_best_observation_time using a dummy Source."""
    src = Source(center.ra, center.dec, 1 * u.Jy)
    return src.get_best_observation_time(telescope=telescope)


# --------------------------------------------------------------------------- #
# simulation
# --------------------------------------------------------------------------- #


def run_simulation(
    ctx: RunContext,
    telescope,
    observation,
    sky_model: SkyModel,
) -> Path:
    """Run InterferometerSimulation and return visibility path."""
    config = ctx.config
    visibility_path = ctx.visibility_path

    if visibility_path.exists():
        if config.overwrite:
            logger.info(f"Overwriting existing {visibility_path}")
            shutil.rmtree(visibility_path)
        else:
            ans = input(f"{visibility_path} exists. Overwrite? (y/n): ")
            if ans.lower() != "y":
                logger.info("User declined overwrite — exiting")
                sys.exit(0)
            shutil.rmtree(visibility_path)

    freq = config.observation.freq_mhz * u.MHz
    fov = compute_fov(config, freq)
    _, _, _, n_channels, delta_freq, _ = build_observation(
        ctx, sky_model.get_center(), telescope
    )

    params = {
        "channel_bandwidth_hz": delta_freq.to(u.Hz).value,
        "station_type": "Gaussian beam",
        "gauss_beam_fwhm_deg": fov.to(u.deg).value,
        "gauss_ref_freq_hz": freq.to(u.Hz).value,
        "use_gpus": False,
    }
    if config.rms:
        params["noise_enable"] = True
        params["noise_freq"] = "Telescope model"
        params["noise_rms"] = "Telescope model"

    simulation = InterferometerSimulation(**params)
    simulation.run_simulation(
        telescope=telescope,
        observation=observation,
        sky=sky_model,
        visibility_path=str(visibility_path),
        backend=SimulatorBackend.OSKAR,
    )
    logger.info(f"Visibilities saved in {visibility_path}")
    ctx.manifest.outputs.append(str(visibility_path.relative_to(ctx.work_dir)))
    return visibility_path


# --------------------------------------------------------------------------- #
# top-level orchestrator
# --------------------------------------------------------------------------- #


def run(config: SimConfig) -> None:
    """Execute the full simulation pipeline from a SimConfig."""
    from .weblog import render_weblog

    t0 = time.time()
    ctx = create_run_context(config)

    try:
        logger.info(f"Telescope : {config.telescope}")
        logger.info(f"Freq      : {config.observation.freq_mhz} MHz")
        logger.info(f"Bandwidth : {config.observation.bandwidth_mhz} MHz")
        logger.info(f"Channels  : {config.observation.n_channels}")
        logger.info(f"Time      : {config.observation.seconds} s")
        logger.info(f"Pixels    : {config.imaging.pixels}")
        logger.info(f"Cleaning  : {config.cleaning}")

        telescope = build_telescope(ctx)
        telescope_png = ctx.work_dir / f"{ctx.work_dir.name}_{config.telescope}_{config.telescope_version or ''}_telescope.png"
        telescope.plot_telescope(file=str(telescope_png))
        ctx.manifest.outputs.append(str(telescope_png.relative_to(ctx.work_dir)))

        freq = config.observation.freq_mhz * u.MHz
        fov = compute_fov(config, freq)
        logger.info(f"FoV       : {fov.to(u.deg).value:.4f} deg")

        sky_model, center = build_sky_model(ctx, fov)
        center = parse_center(config.center, center)
        logger.info(f"Centre    : {center.to_string('hmsdms')}")

        observation, _, bandwidth, n_channels, delta_freq, start_freq = build_observation(ctx, center, telescope)
        logger.info(f"StartFreq : {start_freq.to(u.MHz).value:.3f} MHz")
        logger.info(f"DeltaFreq : {delta_freq.to(u.MHz).value:.3f} MHz")
        logger.info(f"N channels: {n_channels}")

        # phase 1: simulation
        ctx.add_milestone("simulation_started", "started")
        t_phase_a = time.time()
        try:
            visibility_path = run_simulation(ctx, telescope, observation, sky_model)
            ctx.add_milestone("simulation_completed", "completed", elapsed_s=time.time() - t_phase_a)
        except Exception as exc:
            ctx.add_milestone("simulation_failed", "failed", elapsed_s=time.time() - t_phase_a, details=str(exc))
            raise

        # phase 2: imaging
        ctx.add_milestone("imaging_started", "started")
        t_phase_b = time.time()
        try:
            if not config.cleaning:
                run_dirty_imaging(ctx, visibility_path, fov, center)
            else:
                run_wsclean_imaging(ctx, visibility_path, fov)
            ctx.add_milestone("imaging_completed", "completed", elapsed_s=time.time() - t_phase_b, details={"Imager": "dirty" if not config.cleaning else "wsclean"})
        except Exception as exc:
            ctx.add_milestone("imaging_failed", "failed", elapsed_s=time.time() - t_phase_b, details=str(exc))
            raise

        ctx.manifest.mark_completed()
        ctx.save_manifest()

        render_weblog(ctx.manifest, ctx.work_dir)
        logger.info(f"Weblog written to {ctx.weblog_path}")

    except Exception as exc:
        ctx.manifest.mark_failed(str(exc))
        ctx.save_manifest()
        raise

    elapsed = time.time() - t0
    logger.info(f"Done. Elapsed: {elapsed:.1f} s")
