"""Weblog rendering behavior."""

from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.io import fits

from skasim.config import ImgConfig, ObsConfig, SimConfig
from skasim.manifest import RunManifest
from skasim.weblog import _software_versions, render_weblog


def test_weblog_renders_structured_outputs(tmp_path):
    """The weblog displays structured output kinds and paths."""
    manifest = RunManifest(
        run_id="example",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    manifest.add_output("visibility", "visibilities.MS")
    manifest.add_output("log", "example.log")

    html = render_weblog(manifest, tmp_path)

    assert "visibility" in html
    assert "visibilities.MS" in html
    assert "log" in html
    assert "example.log" in html


def test_weblog_skips_missing_image_outputs(tmp_path):
    """Failed runs can render even if an image output record points to a missing file."""
    manifest = RunManifest(
        run_id="missing-image",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    manifest.add_output("image_product", "missing.png", imager="wsclean", role="image")
    manifest.mark_failed("imaging failed before writing image")

    html = render_weblog(manifest, tmp_path)

    assert "imaging failed before writing image" in html
    assert "missing.png" in html


def test_weblog_uses_pipeline_milestone_names_for_durations(tmp_path):
    """Weblog duration cards are derived from actual pipeline milestone names."""
    started = datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc)
    manifest = RunManifest(
        run_id="durations",
        started_at=started,
        config=SimConfig(),
    )
    manifest.add_milestone("simulation_started", "started")
    manifest.milestones[-1].timestamp_utc = started
    manifest.add_milestone("simulation_completed", "completed")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=75)
    manifest.add_milestone("imaging_started", "started")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=80)
    manifest.add_milestone("imaging_completed", "completed")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=90)

    html = render_weblog(manifest, tmp_path)

    assert "simulation 1m 15s" in html
    assert "imaging 10.0s" in html


def test_weblog_groups_wsclean_mfs_products_and_stats(tmp_path):
    """MFS products render as model/clean/residual with peak and RMS stats."""
    manifest = RunManifest(
        run_id="science-products",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    prefix = "science-products_wsclean"
    for role in ("dirty", "model", "image", "residual", "psf"):
        stem = f"{prefix}-MFS-{role if role != 'image' else 'image'}"
        fits_path = tmp_path / f"{stem}.fits"
        png_path = tmp_path / f"{stem}.png"
        data = np.array([[1.0e-3, 2.0e-3], [3.0e-3, 4.0e-3]])
        if role == "residual":
            data = np.array([[0.0, 1.0e-6], [-1.0e-6, 0.0]])
        header = fits.Header()
        header["BMAJ"] = 2.0 / 3600.0
        header["BMIN"] = 1.0 / 3600.0
        header["BPA"] = 35.0
        fits.writeto(fits_path, data, header, overwrite=True)
        png_path.write_bytes(
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
            b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02"
            b"\x00\x00\x00\x90wS\xde\x00\x00\x00\x00IEND\xaeB`\x82"
        )
        preview_role = role
        manifest.add_output(
            "image_product",
            fits_path.name,
            image_product_id=prefix,
            imager="wsclean",
            role=preview_role,
        )
        manifest.add_output(
            "image_product",
            png_path.name,
            image_product_id=prefix,
            imager="wsclean",
            role=f"{preview_role}_preview",
        )

    html = render_weblog(manifest, tmp_path)

    assert "Science Products" in html
    assert "Model" in html
    assert "Clean" in html
    assert "Residual" in html
    assert "Dirty Image" not in html
    assert "View PSF" in html
    assert "Peak:" in html
    assert "4.000 mJy/beam" in html
    assert "RMS:" in html
    assert "BMAJ 2.000 arcsec" in html
    assert "BMIN 1.000 arcsec" in html
    assert "BPA 35.00 deg" in html
    assert "Beam: 2.000 arcsec x 1.000 arcsec" in html
    assert html.count("Beam: 2.000 arcsec x 1.000 arcsec") == 1


def test_weblog_uses_dirty_image_only_when_clean_products_are_absent(tmp_path):
    """OSKAR dirty products render as the fallback science product."""
    manifest = RunManifest(
        run_id="dirty-only",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    fits_path = tmp_path / "dirty-only_dirty.fits"
    png_path = tmp_path / "dirty-only_dirty.png"
    fits.writeto(fits_path, np.ones((2, 2)), overwrite=True)
    png_path.write_bytes(
        b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
        b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02"
        b"\x00\x00\x00\x90wS\xde\x00\x00\x00\x00IEND\xaeB`\x82"
    )
    manifest.add_output(
        "image_product",
        png_path.name,
        image_product_id="dirty-only",
        imager="oskar-dirty",
        role="preview",
    )
    manifest.add_output(
        "image_product",
        fits_path.name,
        image_product_id="dirty-only",
        imager="oskar-dirty",
        role="dirty",
    )

    html = render_weblog(manifest, tmp_path)

    assert "Dirty Image" in html


def test_weblog_renders_observation_imaging_and_cleaning_parameters(tmp_path):
    """Resolved setup tables expose spectral, imaging, and imager parameters."""
    manifest = RunManifest(
        run_id="setup-summary",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(
            observation=ObsConfig(
                frequency_mhz=1300.0,
                bandwidth_mhz=100.0,
                n_channels=8,
                observation_time_s=60,
            ),
            imaging=[ImgConfig(
                imager="wsclean",
                pixels=1024,
                fov_deg=1.0,
                robust=-0.5,
                wsclean_command="wsclean",
                clean_iterations=500,
            )],
        ),
    )
    manifest.add_milestone(
        "observation_configured",
        "completed",
        details={
            "min_frequency_mhz": 1250.0,
            "max_frequency_mhz": 1350.0,
            "n_timesteps": 7,
            "phase_center_ra_deg": 150.0,
            "phase_center_dec_deg": 2.0,
        },
    )
    manifest.add_output("visibility", "visibilities.MS")
    manifest.add_output(
        "image_product",
        "setup-summary_wsclean-MFS-image.fits",
        image_product_id="setup-summary_wsclean",
        imager="wsclean",
        role="image",
    )

    html = render_weblog(manifest, tmp_path)

    assert "Observation & Telescope" in html
    assert "Frequency Setup" not in html
    assert "Frequency Range" in html
    assert "1250-1350 MHz" in html
    assert "Central Frequency" in html
    assert "1300 MHz" in html
    assert "Channels" in html
    assert "8 x 12.5 MHz" in html
    assert "12.5 MHz" in html
    assert "Total Bandwidth" in html
    assert "100 MHz" in html
    assert "Imaging Setup" in html
    assert "1024 x 1024 pixels" in html
    assert "Pixels</th>" not in html
    assert "Total FoV" in html
    assert "1 deg" in html
    assert "Pixel Size" in html
    assert "3.5156 arcsec" in html
    assert "Cleaning & Imager Parameters" in html
    assert "Clean iterations" in html
    assert "500" in html
    assert "Major-cycle gain" in html
    assert "0.8" in html
    assert "Auto threshold" in html
    assert "Auto mask" in html
    assert "Channels out" in html
    assert "Output prefix" in html
    assert "setup-summary_wsclean" in html
    assert "Visibility input" in html
    assert "visibilities.MS" in html


def test_weblog_renders_antenna_count_in_telescope_section(tmp_path):
    """Telescope metadata exposes the station/antenna count with a clear label."""
    manifest = RunManifest(
        run_id="antennas",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(telescope="MeerKAT"),
    )
    manifest.add_milestone(
        "telescope_built",
        "completed",
        details={"n_stations": 64},
    )

    html = render_weblog(manifest, tmp_path)

    assert "Antennas" in html
    assert "64" in html


def test_weblog_uses_known_antenna_count_for_older_meerkat_manifests(tmp_path):
    """Older manifests without telescope-built counts still show known MeerKAT antennas."""
    manifest = RunManifest(
        run_id="old-meerkat",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(telescope="MeerKAT"),
    )
    manifest.add_milestone(
        "telescope_built",
        "completed",
        details={"name": "MeerKAT", "version": None},
    )

    html = render_weblog(manifest, tmp_path)

    assert "Antennas" in html
    assert "64" in html


def test_weblog_renders_software_versions(tmp_path):
    """Weblog includes a concise reproducibility software version table."""
    manifest = RunManifest(
        run_id="versions",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )

    html = render_weblog(manifest, tmp_path)

    assert "Software Versions" in html
    for name, _ in _software_versions():
        assert name in html
