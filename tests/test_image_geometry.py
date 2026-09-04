"""Tests for the backend-independent image-geometry contract."""

import math
import warnings

import astropy.units as u
import pytest
from pydantic import ValidationError

from skasim import pipeline
from skasim.config import ImgConfig, ObsConfig, SimConfig
from skasim.image_geometry import (
    GEOMETRY_RELATIVE_TOLERANCE,
    LEGACY_IMAGE_PIXELS,
    MAX_IMAGE_PIXELS,
    resolve_image_geometry,
)
from skasim.manifest import create_run_context
from skasim.weblog import render_weblog


def _resolve(**requested):
    return resolve_image_geometry(
        fov_deg=requested.get("fov_deg"),
        pixels=requested.get("pixels"),
        cell_size_arcsec=requested.get("cell_size_arcsec"),
        diffraction_fov_deg=1.0,
        theoretical_beam_arcsec=10.0,
        reference_frequency_hz=1.5e9,
    )


@pytest.mark.parametrize(
    "requested, expected_pixels, expected_fov, expected_cell, legacy",
    [
        ({}, LEGACY_IMAGE_PIXELS, 1.0, 3600.0 / 512, True),
        ({"fov_deg": 1.0}, 2048, 1.0, 3600.0 / 2048, False),
        ({"pixels": 900}, 900, 1.0, 4.0, False),
        ({"cell_size_arcsec": 4.0}, 1024, 1024 * 4.0 / 3600.0, 4.0, False),
        ({"fov_deg": 1.0, "pixels": 1000}, 1000, 1.0, 3.6, False),
        (
            {"fov_deg": 1.0, "cell_size_arcsec": 4.0},
            1024,
            1024 * 4.0 / 3600.0,
            4.0,
            False,
        ),
        (
            {"pixels": 900, "cell_size_arcsec": 4.0},
            900,
            1.0,
            4.0,
            False,
        ),
        (
            {"fov_deg": 1.0, "pixels": 900, "cell_size_arcsec": 4.0},
            900,
            1.0,
            4.0,
            False,
        ),
    ],
)
def test_resolves_every_parameter_combination(
    requested, expected_pixels, expected_fov, expected_cell, legacy
):
    geometry = _resolve(**requested)

    assert geometry.effective_pixels == expected_pixels
    assert geometry.effective_fov_deg == pytest.approx(expected_fov)
    assert geometry.effective_cell_size_arcsec == pytest.approx(expected_cell)
    assert geometry.legacy_fallback is legacy


def test_fov_only_targets_five_pixels_per_beam_subject_to_backend_minimum():
    geometry = _resolve(fov_deg=(10.0 / 5.0) / 3600.0 * 7)

    assert geometry.effective_pixels == 64
    assert geometry.effective_fov_deg == pytest.approx((10.0 / 5.0) / 3600.0 * 7)
    assert geometry.pixels_rounded_up is True
    assert geometry.pixels_per_beam == pytest.approx(64 / 7 * 5)


def test_fov_only_respects_maximum_dimension():
    huge_fov = MAX_IMAGE_PIXELS * 2.0 / 5.0 * 10.0 / 3600.0
    geometry = _resolve(fov_deg=huge_fov)

    assert geometry.effective_pixels == MAX_IMAGE_PIXELS
    assert "exceeding the maximum" in " ".join(geometry.warnings)
    assert geometry.effective_fov_deg == pytest.approx(huge_fov)
    assert geometry.pixels_rounded_up is False


def test_cell_size_only_preserves_cell_and_does_not_crop_diffraction_fov():
    geometry = _resolve(cell_size_arcsec=4.0)
    assert geometry.effective_pixels == 1024
    assert geometry.effective_cell_size_arcsec == pytest.approx(4.0)
    assert geometry.effective_fov_deg >= 1.0


def test_fov_and_cell_size_preserves_cell_and_does_not_crop_requested_fov():
    geometry = _resolve(fov_deg=0.5, cell_size_arcsec=2.0)
    assert geometry.effective_pixels == 1024
    assert geometry.effective_cell_size_arcsec == pytest.approx(2.0)
    assert geometry.effective_fov_deg >= 0.5


def test_explicit_non_power_of_two_is_preserved():
    assert _resolve(fov_deg=1.0, pixels=777).effective_pixels == 777


@pytest.mark.parametrize(
    "fov_deg",
    [1.0 / (1.0 + GEOMETRY_RELATIVE_TOLERANCE), 1.0],
)
def test_three_value_consistency_accepts_tolerance_boundary(fov_deg):
    geometry = _resolve(fov_deg=fov_deg, pixels=900, cell_size_arcsec=4.0)
    assert geometry.effective_pixels == 900


def test_three_value_consistency_rejects_beyond_tolerance():
    with pytest.raises(ValueError, match="inconsistent image geometry"):
        _resolve(
            fov_deg=1.0 - GEOMETRY_RELATIVE_TOLERANCE * 1.01,
            pixels=900,
            cell_size_arcsec=4.0,
        )


def test_config_rejects_inconsistent_fully_specified_geometry():
    with pytest.raises(ValidationError, match="inconsistent image geometry"):
        ImgConfig(fov_deg=1.0, pixels=900, cell_size_arcsec=4.1)


def test_explicit_default_pixel_count_remains_authoritative():
    config = ImgConfig(fov_deg=1.0, pixels=512)

    geometry = config.resolve_geometry(1.0, 10.0, 1.5e9)

    assert geometry.requested_pixels == 512
    assert geometry.effective_pixels == 512


def test_explicit_default_pixel_count_does_not_bypass_triplet_validation():
    with pytest.raises(ValidationError, match="inconsistent image geometry"):
        ImgConfig(fov_deg=1.0, pixels=512, cell_size_arcsec=99.0)


def test_explicit_maximum_dimension_is_allowed_and_excess_is_rejected():
    assert _resolve(pixels=16_384).effective_pixels == 16_384
    with pytest.raises(ValidationError, match="pixels must be <= 16384"):
        ImgConfig(pixels=16_385)


def test_derived_dimension_at_limit_and_above_limit():
    at_limit = _resolve(fov_deg=16_384 / 3600, cell_size_arcsec=1.0)
    assert at_limit.effective_pixels == 16_384
    with pytest.raises(ValueError, match="exceeds the maximum"):
        _resolve(fov_deg=16_385 / 3600, cell_size_arcsec=1.0)


@pytest.mark.parametrize(
    "cell_size, phrase",
    [(5.0, "under-samples"), (0.5, "over-samples")],
)
def test_sampling_warnings(cell_size, phrase):
    geometry = _resolve(pixels=100, cell_size_arcsec=cell_size)
    assert any(phrase in message for message in geometry.warnings)


def test_legacy_config_emits_validation_deprecation_and_resolver_warning():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ImgConfig()
    assert any(item.category is DeprecationWarning for item in caught)
    assert any("Legacy image geometry" in message for message in _resolve().warnings)


def test_result_metadata_is_serializable_and_complete():
    record = _resolve(fov_deg=1.0).as_dict()
    assert math.isfinite(record["theoretical_beam_arcsec"])
    assert record["reference_frequency_hz"] == 1.5e9
    assert "requested_fov_deg" in record
    assert "effective_cell_size_arcsec" in record


def test_pipeline_resolves_blocks_independently_and_records_weblog(tmp_path):
    class TelescopeFixture:
        def max_baseline(self):
            return 10_000.0

    config = SimConfig(
        output_dir=str(tmp_path / "geometry-run"),
        observation=ObsConfig(
            frequency_mhz=1000.0,
            bandwidth_mhz=200.0,
            n_channels=2,
        ),
        imaging=[
            ImgConfig(tag="wide", fov_deg=1.0),
            ImgConfig(tag="compact", pixels=1000, cell_size_arcsec=10.0),
        ],
    )
    ctx = create_run_context(config)

    _, _, _, geometries = pipeline._resolve_run_geometry(ctx, TelescopeFixture())

    assert geometries["wide"].effective_pixels != geometries["compact"].effective_pixels
    assert geometries["wide"].theoretical_beam_arcsec == pytest.approx(
        geometries["compact"].theoretical_beam_arcsec
    )
    assert geometries["wide"].reference_frequency_hz == pytest.approx(1.0e9)
    expected_beam = (1.0e9 * u.Hz).to(u.m, equivalencies=u.spectral()).value / 10_000
    expected_beam = (expected_beam * u.rad).to(u.arcsec).value
    assert geometries["wide"].theoretical_beam_arcsec == pytest.approx(expected_beam)
    geometry_milestone = next(
        item
        for item in ctx.manifest.milestones
        if item.name == "image_geometry_resolved"
    )
    assert geometry_milestone.details["simulation_fov_deg"] == pytest.approx(
        geometries["compact"].effective_fov_deg
    )
    assert (
        geometry_milestone.details["diffraction_fov_frequency_label"] == "band centre"
    )
    central_diffraction_fov = (
        pipeline.compute_fov(config.telescope, None, 1000.0 * u.MHz).to(u.deg).value
    )
    assert geometry_milestone.details["diffraction_fov_deg"] == pytest.approx(
        central_diffraction_fov
    )

    html = render_weblog(ctx.manifest, ctx.work_dir)
    assert "4096 x 4096" in html
    assert "Requested Geometry" in html
    assert "Theoretical Beam" in html
    assert "1000 MHz" in html
    assert "band centre" in html
    assert "Pixels per Beam" in html
    assert "Power-of-two Rounding" in html
    assert "Legacy Fallback" in html
    assert "Geometry Warnings" in html
