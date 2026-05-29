"""pipeline.py — end-to-end interferometric simulation orchestrator."""

from __future__ import annotations

import json
import os
import pickle
import shutil
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from loguru import logger

from .config import SimConfig
from .imaging import run_dirty_imaging, run_wsclean_imaging, write_uv_coverage_plot
from .manifest import RunContext, create_run_context
from .fits_helper import FitsCatalogLoader
from .runtime import require_karabo_module
from .sky import SkyModel, Source
from .utils import get_diameter, mapping_unit

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
    simulator_backend = require_karabo_module("karabo.simulator_backend")
    telescope_module = require_karabo_module("karabo.simulation.telescope")
    kwargs: dict = {"backend": simulator_backend.SimulatorBackend.OSKAR}
    if config.telescope_version is not None:
        kwargs["version"] = resolve_telescope_version(
            telescope_module,
            config.telescope,
            config.telescope_version,
        )
        logger.info(f"Telescope {config.telescope}  version={config.telescope_version}")
    else:
        logger.info(f"Telescope {config.telescope}  (no version)")
    telescope = telescope_module.Telescope.constructor(config.telescope, **kwargs)
    details = {"name": config.telescope, "version": config.telescope_version}
    try:
        if hasattr(telescope, "antennas") and telescope.antennas is not None:
            details["n_stations"] = len(telescope.antennas)
        elif hasattr(telescope, "num_antennas"):
            details["n_stations"] = telescope.num_antennas
        elif hasattr(telescope, "num_stations"):
            details["n_stations"] = telescope.num_stations
    except Exception:
        pass
    ctx.add_milestone(
        "telescope_built",
        "completed",
        details=details,
    )
    return telescope


def resolve_telescope_version(telescope_module, telescope: str, version: str):
    """Resolve a CLI/config telescope version string to Karabo's enum member."""
    version_enum = getattr(telescope_module, "OSKAR_TELESCOPE_TO_VERSIONS", {}).get(
        telescope
    )
    if version_enum is None:
        return version
    if not isinstance(version, str):
        return version
    for candidate in version_enum:
        if version in (candidate.name, candidate.value):
            return candidate
    accepted = ", ".join(candidate.name for candidate in version_enum)
    raise ValueError(
        f"Unsupported version {version!r} for telescope {telescope!r}. "
        f"Accepted versions: {accepted}"
    )


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
    flux_scale: float = 1.0,
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

    # json catalog
    if ext == ".json":
        with open(fpath, "r") as fh:
            data = json.load(fh)
        sources: List[Source] = []
        for item in data:
            src = Source.from_json(item)
            if flux_scale != 1.0:
                src.I *= flux_scale
            if src.ref_freq == 0:
                hz = frequency.to(u.Hz).value if frequency is not None else 0
                src.ref_freq = hz * u.Hz
            sources.append(src)
        if not sources:
            raise ValueError(f"No sources found in JSON {fpath}")
        sky_model = SkyModel()
        arr = np.array([s.to_sky_model(reduced_form=False) for s in sources])
        sky_model.add_point_sources(arr)
        sky_model.get_center()
        logger.info(f"Loaded {len(sources)} sources from JSON {fpath}")
        return sky_model

    # fits table or image
    if ext in (".fits", ".fit"):
        return _load_sky_from_fits(
            fpath, column_mapping, flux_scale, frequency
        )

    raise ValueError(f"Unsupported sky-file extension: {ext}")



def _load_sky_from_fits(
    fpath: str,
    column_mapping: str,
    flux_scale: float,
    frequency: Optional[u.Quantity],
) -> SkyModel:
    """try Karabo's get_sky_model_from_fits; fallback to our own loader when columns lack TUNIT."""
    loader = FitsCatalogLoader(
        fpath=fpath,
        column_mapping=column_mapping,
        scale_I=flux_scale,
        ref_freq_hz=None,
        frequency=frequency,
    )
    return loader.load()


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
            flux_scale=config.flux_scale,
            frequency=config.observation.frequency_mhz * u.MHz,
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
        ctx.manifest.add_output(
            "sky_model",
            str(ctx.sky_file_resolved),
            metadata={"format": config.sky_format, "n_sources": n_srcs},
        )
        return sky_model, center

    # 2) built-in catalog
    if config.catalog is not None:
        if config.catalog == "MIGHTEE":
            logger.info("Loading MIGHTEE catalog")
            if not hasattr(SkyModel, "get_MIGHTEE_Sky"):
                require_karabo_module("karabo.simulation.sky_model")
            sky_model = SkyModel.get_MIGHTEE_Sky()
            fmt = "MIGHTEE"
        elif config.catalog == "GLEAM":
            logger.info("Loading GLEAM catalog")
            if not hasattr(SkyModel, "get_GLEAM_Sky"):
                require_karabo_module("karabo.simulation.sky_model")
            sky_model = SkyModel.get_GLEAM_Sky()
            fmt = "GLEAM"
        elif config.catalog == "SKAMid":
            skamid_path = Path("SKAMid_B1_8h_v3.fits").resolve()
            if skamid_path.exists():
                logger.info(f"Loading SKAMid catalog {skamid_path}")
                sky_model = SkyModel.get_sky_model_from_fits(fits_file=str(skamid_path))
                fmt = "SKAMid"
            else:
                logger.info(f"SKAMid catalog not found at {skamid_path}")
                raise FileNotFoundError(str(skamid_path))
        else:
            raise ValueError(f"Catalog {config.catalog} not available")
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
    intensities = [i * u.Jy for i in config.source_flux_jy]
    stokes_q = config.stokes_q_jy or [0.0] * len(intensities)
    stokes_u = config.stokes_u_jy or [0.0] * len(intensities)
    stokes_v = config.stokes_v_jy or [0.0] * len(intensities)
    n_sources = len(intensities)
    sources: List[Source] = []
    for idx in range(n_sources):
        if idx == 0:
            src = source_ref
            src.I = intensities[idx]
            src.Q = stokes_q[idx] * u.Jy
            src.U = stokes_u[idx] * u.Jy
            src.V = stokes_v[idx] * u.Jy
        else:
            x_coord = np.random.uniform(-fov.value / 2, fov.value / 2) * 0.8 * u.rad
            y_coord = np.random.uniform(-fov.value / 2, fov.value / 2) * 0.8 * u.rad
            src = Source(
                source_ref.ra + x_coord,
                source_ref.dec + y_coord,
                intensities[idx],
                Q=stokes_q[idx] * u.Jy,
                U=stokes_u[idx] * u.Jy,
                V=stokes_v[idx] * u.Jy,
            )
        sources.append(src)

    sky_model = SkyModel()
    has_polarization = any(
        value != 0.0 for values in (stokes_q, stokes_u, stokes_v) for value in values
    )
    arr = np.array(
        [s.to_sky_model(reduced_form=not has_polarization) for s in sources]
    )
    sky_model.add_point_sources(arr)
    center = sky_model.get_center()
    ctx.add_milestone(
        "sky_model_loaded",
        "completed",
        details={"format": "random", "n_sources": n_sources, "reference": "HCG16"},
    )
    return sky_model, center


def build_observation(
    ctx: RunContext,
    center: SkyCoord,
    telescope,
) -> tuple:
    """Return (observation, frequency, bandwidth, n_channels, delta_freq, start_freq)."""
    observation_module = require_karabo_module("karabo.simulation.observation")
    config = ctx.config
    obs = config.observation
    freq = obs.frequency_mhz * u.MHz
    bw_mhz = obs.bandwidth_mhz
    n_channels = obs.n_channels
    df_mhz = obs.channel_width_mhz
    bandwidth = bw_mhz * u.MHz
    delta_freq = df_mhz * u.MHz
    seconds = config.observation.observation_time_s

    start_freq = freq - n_channels * delta_freq / 2
    end_freq = start_freq + n_channels * delta_freq

    # best observation time (culmination)
    obs_time = source_ref_get_best_observation_time(center, telescope)
    n_timesteps = max(1, int(seconds / 7.997))

    observation = observation_module.Observation(
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
            "frequency_mhz": obs.frequency_mhz,
            "min_frequency_mhz": start_freq.to(u.MHz).value,
            "max_frequency_mhz": end_freq.to(u.MHz).value,
            "bandwidth_mhz": bw_mhz,
            "n_channels": n_channels,
            "channel_width_mhz": df_mhz,
            "observation_time_s": seconds,
            "n_timesteps": n_timesteps,
            "phase_center_ra_deg": center.ra.to(u.deg).value,
            "phase_center_dec_deg": center.dec.to(u.deg).value,
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
    interferometer_module = require_karabo_module("karabo.simulation.interferometer")
    simulator_backend = require_karabo_module("karabo.simulator_backend")
    config = ctx.config
    visibility_path = ctx.visibility_path

    if visibility_path.exists():
        if config.overwrite:
            logger.info(f"Overwriting existing {visibility_path}")
            shutil.rmtree(visibility_path)
        else:
            raise FileExistsError(
                f"{visibility_path} already exists. Use --overwrite to replace it."
            )

    freq = config.observation.frequency_mhz * u.MHz
    fov = compute_fov(config, freq)
    delta_freq = config.observation.channel_width_mhz * u.MHz

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

    simulation = interferometer_module.InterferometerSimulation(**params)
    simulation.run_simulation(
        telescope=telescope,
        observation=observation,
        sky=sky_model,
        visibility_path=str(visibility_path),
        backend=simulator_backend.SimulatorBackend.OSKAR,
    )
    logger.info(f"Visibilities saved in {visibility_path}")
    ctx.manifest.add_output(
        "visibility",
        str(visibility_path.relative_to(ctx.work_dir)),
    )
    return visibility_path


# --------------------------------------------------------------------------- #
# top-level orchestrator
# --------------------------------------------------------------------------- #


def run(config: SimConfig) -> None:
    """Execute the full simulation pipeline from a SimConfig."""
    from .weblog import render_weblog

    import matplotlib

    matplotlib.use("Agg", force=True)

    t0 = time.time()
    ctx = create_run_context(config)

    try:
        logger.info(f"Telescope : {config.telescope}")
        logger.info(f"Frequency : {config.observation.frequency_mhz} MHz")
        logger.info(f"Bandwidth : {config.observation.bandwidth_mhz} MHz")
        logger.info(f"Channels  : {config.observation.n_channels}")
        logger.info(f"Obs time  : {config.observation.observation_time_s} s")
        logger.info(f"Pixels    : {config.imaging.pixels}")
        logger.info(f"Imager    : {config.imaging.imager}")

        telescope = build_telescope(ctx)
        telescope_png = ctx.work_dir / f"{ctx.work_dir.name}_{config.telescope}_{config.telescope_version or ''}_telescope.png"
        try:
            telescope.plot_telescope(file=str(telescope_png))
        finally:
            import matplotlib.pyplot as plt

            plt.close("all")
        ctx.manifest.add_output(
            "plot",
            str(telescope_png.relative_to(ctx.work_dir)),
            role="telescope",
        )

        freq = config.observation.frequency_mhz * u.MHz
        fov = compute_fov(config, freq)
        logger.info(f"FoV       : {fov.to(u.deg).value:.4f} deg")

        sky_model, center = build_sky_model(ctx, fov)
        center = parse_center(config.center, center)
        logger.info(f"Centre    : {center.to_string('hmsdms')}")
        try:
            from .imaging import write_sky_model_previews

            for path, role in write_sky_model_previews(
                sky_model,
                center,
                fov,
                ctx.work_dir,
                ctx.work_dir.name,
            ):
                ctx.manifest.add_output("plot", path, role=role)
        except Exception as exc:
            logger.warning(f"Sky model previews failed: {show_exc(exc)}")

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
            ctx.add_milestone("simulation_failed", "failed", elapsed_s=time.time() - t_phase_a, details={"error": str(exc)})
            raise

        if config.imaging.uv_coverage:
            ctx.add_milestone("uv_coverage_started", "started")
            t_phase_uv = time.time()
            try:
                uv_coverage_path = write_uv_coverage_plot(ctx, visibility_path)
                ctx.add_milestone(
                    "uv_coverage_completed",
                    "completed",
                    elapsed_s=time.time() - t_phase_uv,
                    details={"path": str(uv_coverage_path.relative_to(ctx.work_dir))},
                )
            except Exception as exc:
                ctx.add_milestone(
                    "uv_coverage_failed",
                    "failed",
                    elapsed_s=time.time() - t_phase_uv,
                    details={"error": str(exc)},
                )
                raise

        # phase 2: imaging
        ctx.add_milestone("imaging_started", "started")
        t_phase_b = time.time()
        try:
            if config.imaging.imager == "oskar-dirty":
                run_dirty_imaging(ctx, visibility_path, fov, center)
            else:
                run_wsclean_imaging(ctx, visibility_path, fov)
            ctx.add_milestone(
                "imaging_completed",
                "completed",
                elapsed_s=time.time() - t_phase_b,
                details={"Imager": config.imaging.imager},
            )
        except Exception as exc:
            ctx.add_milestone("imaging_failed", "failed", elapsed_s=time.time() - t_phase_b, details={"error": str(exc)})
            raise

        ctx.manifest.mark_completed()
        ctx.manifest.add_output("weblog", ctx.weblog_path.name)
        ctx.save_manifest()

        render_weblog(ctx.manifest, ctx.work_dir)
        logger.info(f"Weblog written to {ctx.weblog_path}")

    except Exception as exc:
        ctx.manifest.mark_failed(str(exc))
        ctx.manifest.add_output("weblog", ctx.weblog_path.name)
        ctx.save_manifest()
        render_weblog(ctx.manifest, ctx.work_dir)
        logger.info(f"Failure weblog written to {ctx.weblog_path}")
        raise

    elapsed = time.time() - t0
    logger.info(f"Done. Elapsed: {elapsed:.1f} s")
