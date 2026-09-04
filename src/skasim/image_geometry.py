"""Pure image-geometry resolution independent of imaging backends."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Optional, Tuple

MAX_IMAGE_PIXELS = 16_384
MIN_IMAGE_PIXELS = 64
GEOMETRY_RELATIVE_TOLERANCE = 1e-6
LEGACY_IMAGE_PIXELS = 512
TARGET_PIXELS_PER_BEAM = 5.0


@dataclass(frozen=True)
class ImageGeometry:
    """Requested and effective geometry for one imaging block."""

    requested_fov_deg: Optional[float]
    requested_pixels: Optional[int]
    requested_cell_size_arcsec: Optional[float]
    effective_fov_deg: float
    effective_pixels: int
    effective_cell_size_arcsec: float
    theoretical_beam_arcsec: float
    reference_frequency_hz: float
    pixels_per_beam: float
    pixels_rounded_up: bool
    legacy_fallback: bool
    warnings: Tuple[str, ...]

    def as_dict(self) -> dict:
        """Return a JSON-serializable run-record representation."""
        value = asdict(self)
        value["warnings"] = list(self.warnings)
        return value


def validate_geometry_triplet(
    fov_deg: float, pixels: int, cell_size_arcsec: float
) -> None:
    """Validate a fully specified triplet using the public tolerance contract."""
    implied_fov_deg = pixels * cell_size_arcsec / 3600.0
    if not math.isclose(
        fov_deg,
        implied_fov_deg,
        rel_tol=GEOMETRY_RELATIVE_TOLERANCE,
        abs_tol=0.0,
    ):
        raise ValueError(
            "inconsistent image geometry: fov_deg must equal "
            "pixels * cell_size_arcsec / 3600 within a relative tolerance of 1e-6"
        )


def resolve_image_geometry(
    *,
    fov_deg: Optional[float],
    pixels: Optional[int],
    cell_size_arcsec: Optional[float],
    diffraction_fov_deg: float,
    theoretical_beam_arcsec: float,
    reference_frequency_hz: float,
) -> ImageGeometry:
    """Resolve one geometry triplet and report its sampling characteristics.

    Automatically selected dimensions are rounded upward to powers of two.
    Explicit cell sizes remain fixed; their effective field may expand so the
    requested or diffraction-limited field is not cropped. For field-only
    requests, the field remains fixed and the effective cell size is adjusted.
    """
    _require_positive("diffraction_fov_deg", diffraction_fov_deg)
    _require_positive("theoretical_beam_arcsec", theoretical_beam_arcsec)
    _require_positive("reference_frequency_hz", reference_frequency_hz)
    if fov_deg is not None:
        _require_positive("fov_deg", fov_deg)
    if cell_size_arcsec is not None:
        _require_positive("cell_size_arcsec", cell_size_arcsec)
    if pixels is not None:
        if pixels < 1:
            raise ValueError("pixels must be >= 1")
        _check_dimension(pixels)

    requested = (fov_deg, pixels, cell_size_arcsec)
    defined = sum(value is not None for value in requested)
    legacy = defined == 0
    rounded = False
    warning_messages: list[str] = []

    if defined == 3:
        assert fov_deg is not None and pixels is not None
        assert cell_size_arcsec is not None
        validate_geometry_triplet(fov_deg, pixels, cell_size_arcsec)
        effective_fov = fov_deg
        effective_pixels = pixels
        effective_cell = cell_size_arcsec
    elif defined == 0:
        effective_fov = diffraction_fov_deg
        effective_pixels = LEGACY_IMAGE_PIXELS
        effective_cell = effective_fov * 3600.0 / effective_pixels
        warning_messages.append(
            "Legacy image geometry is in use; explicitly set fov_deg, pixels, "
            "or cell_size_arcsec."
        )
    elif pixels is not None and cell_size_arcsec is not None:
        effective_pixels = pixels
        effective_cell = cell_size_arcsec
        effective_fov = effective_pixels * effective_cell / 3600.0
    elif fov_deg is not None and pixels is not None:
        effective_fov = fov_deg
        effective_pixels = pixels
        effective_cell = effective_fov * 3600.0 / effective_pixels
    elif fov_deg is not None and cell_size_arcsec is not None:
        # Preserve the explicit cell size and round the derived dimension up.
        # The resulting field may expand, but never crops the requested field.
        effective_pixels, rounded = _derive_pixels(fov_deg * 3600.0 / cell_size_arcsec)
        effective_cell = cell_size_arcsec
        effective_fov = effective_pixels * effective_cell / 3600.0
    elif fov_deg is not None:
        target_cell = theoretical_beam_arcsec / TARGET_PIXELS_PER_BEAM
        required_pixels = fov_deg * 3600.0 / target_cell
        if required_pixels > MAX_IMAGE_PIXELS:
            warning_messages.append(
                f"Requested fov_deg={fov_deg} would require {int(math.ceil(required_pixels))} "
                f"pixels at {TARGET_PIXELS_PER_BEAM} pixels per beam, exceeding the maximum "
                f"of {MAX_IMAGE_PIXELS}. Falling back to the maximum dimension; the effective "
                f"cell size will be smaller and may over-sample the beam."
            )
            effective_pixels = MAX_IMAGE_PIXELS
            rounded = False
        else:
            effective_pixels, rounded = _derive_pixels(required_pixels)
        effective_fov = fov_deg
        effective_cell = effective_fov * 3600.0 / effective_pixels
    elif pixels is not None:
        effective_fov = diffraction_fov_deg
        effective_pixels = pixels
        effective_cell = effective_fov * 3600.0 / effective_pixels
    else:
        # Preserve the explicit cell size. Rounding the derived dimension up
        # may expand the diffraction-limited fallback field, but never crops it.
        assert cell_size_arcsec is not None
        effective_pixels, rounded = _derive_pixels(
            diffraction_fov_deg * 3600.0 / cell_size_arcsec
        )
        effective_cell = cell_size_arcsec
        effective_fov = effective_pixels * effective_cell / 3600.0

    _check_dimension(effective_pixels)
    pixels_per_beam = theoretical_beam_arcsec / effective_cell
    if pixels_per_beam < 3.0:
        warning_messages.append(
            f"Image geometry under-samples the theoretical beam "
            f"({pixels_per_beam:.3g} pixels per beam; recommended minimum is 3)."
        )
    elif pixels_per_beam > 10.0:
        warning_messages.append(
            f"Image geometry over-samples the theoretical beam "
            f"({pixels_per_beam:.3g} pixels per beam; recommended maximum is 10)."
        )

    return ImageGeometry(
        requested_fov_deg=float(fov_deg) if fov_deg is not None else None,
        requested_pixels=int(pixels) if pixels is not None else None,
        requested_cell_size_arcsec=(
            float(cell_size_arcsec) if cell_size_arcsec is not None else None
        ),
        effective_fov_deg=float(effective_fov),
        effective_pixels=int(effective_pixels),
        effective_cell_size_arcsec=float(effective_cell),
        theoretical_beam_arcsec=float(theoretical_beam_arcsec),
        reference_frequency_hz=float(reference_frequency_hz),
        pixels_per_beam=float(pixels_per_beam),
        pixels_rounded_up=bool(rounded),
        legacy_fallback=legacy,
        warnings=tuple(warning_messages),
    )


def _derive_pixels(required_pixels: float) -> Tuple[int, bool]:
    required = max(MIN_IMAGE_PIXELS, math.ceil(required_pixels))
    derived = 1 << (required - 1).bit_length()
    _check_dimension(derived)
    return derived, bool(derived != required_pixels)


def _check_dimension(pixels: int) -> None:
    if pixels > MAX_IMAGE_PIXELS:
        raise ValueError(
            f"image dimension {pixels} exceeds the maximum of {MAX_IMAGE_PIXELS} pixels"
        )


def _require_positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite value greater than zero")
