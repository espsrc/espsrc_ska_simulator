"""Image-model validation, preview, and CASA injection helpers."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from loguru import logger

from .config import (
    ComponentSkyModelEntry,
    ContinuumIAlphaModelEntry,
    ModelEntry,
    SimConfig,
    StaticStokesMapsModelEntry,
)
from .imaging import write_fits_preview
from .manifest import RunContext


@dataclass(frozen=True)
class FitsImageInfo:
    """Small summary of one accepted FITS image model plane."""

    path: Path
    spatial_shape: tuple[int, int]
    unit: str | None
    celestial_header: dict[str, object]
    center: SkyCoord | None


@dataclass(frozen=True)
class CasaModelProduct:
    """CASA-ready model product generated from one model entry."""

    model_paths: list[Path]
    nterms: int
    reffreq: str
    intermediates: list[Path]


def component_model_entries(config: SimConfig) -> list[ComponentSkyModelEntry]:
    return [
        entry
        for entry in config.models
        if isinstance(entry, ComponentSkyModelEntry)
    ]


def image_model_entries(config: SimConfig) -> list[ModelEntry]:
    return [
        entry
        for entry in config.models
        if isinstance(entry, (ContinuumIAlphaModelEntry, StaticStokesMapsModelEntry))
    ]


def image_model_center(entries: list[ModelEntry]) -> SkyCoord | None:
    """Return the centre of the first image model with usable celestial WCS."""
    for entry in entries:
        path = primary_model_fits_path(entry)
        if path is None:
            continue
        try:
            info = read_fits_image_info(path)
        except Exception:
            continue
        if info.center is not None:
            return info.center
    return None


def primary_model_fits_path(entry: ModelEntry) -> Path | None:
    """Return the representative FITS image for previews and phase-centre inference."""
    if isinstance(entry, ContinuumIAlphaModelEntry):
        return Path(entry.stokes_i).expanduser().resolve()
    if isinstance(entry, StaticStokesMapsModelEntry):
        for value in (entry.stokes_i, entry.stokes_q, entry.stokes_u, entry.stokes_v):
            if value:
                return Path(value).expanduser().resolve()
    return None


def read_fits_image_info(path: Path) -> FitsImageInfo:
    """Read FITS image metadata used by validation and reporting."""
    with fits.open(path) as hdul:
        hdu = hdul[0]
        if hdu.data is None:
            raise ValueError(f"{path} has no image data")
        data = np.asarray(hdu.data).squeeze()
        if data.ndim < 2:
            raise ValueError(f"{path} is not a spatial FITS image")
        spatial_shape = tuple(int(v) for v in data.shape[-2:])
        header = hdu.header.copy()
        unit = header.get("BUNIT")

    try:
        wcs = WCS(header).celestial
        celestial_header = dict(wcs.to_header())
        center_y = (spatial_shape[0] - 1) / 2.0
        center_x = (spatial_shape[1] - 1) / 2.0
        center = wcs.pixel_to_world(center_x, center_y)
        if not isinstance(center, SkyCoord):
            center = None
    except Exception:
        celestial_header = {}
        center = None

    return FitsImageInfo(
        path=path,
        spatial_shape=spatial_shape,
        unit=unit,
        celestial_header=celestial_header,
        center=center,
    )


def validate_continuum_i_alpha(entry: ContinuumIAlphaModelEntry) -> dict:
    """Validate the continuum image contract and return report metadata."""
    stokes_info = read_fits_image_info(Path(entry.stokes_i).expanduser().resolve())
    alpha_info = read_fits_image_info(Path(entry.alpha).expanduser().resolve())
    if stokes_info.spatial_shape != alpha_info.spatial_shape:
        raise ValueError(
            "continuum_i_alpha requires matching spatial dimensions: "
            f"{stokes_info.path} has {stokes_info.spatial_shape}, "
            f"{alpha_info.path} has {alpha_info.spatial_shape}"
        )
    if stokes_info.celestial_header != alpha_info.celestial_header:
        raise ValueError("continuum_i_alpha requires matching celestial WCS.")
    unit = (stokes_info.unit or "").strip().lower()
    if unit not in {"jy/pixel", "jy pix-1", "jy/pix", "jy"}:
        raise ValueError(
            f"{stokes_info.path} must declare Jy/pixel-compatible BUNIT; "
            f"found {stokes_info.unit!r}"
        )
    alpha_unit = (alpha_info.unit or "").strip().lower()
    if alpha_unit not in {"", "1", "dimensionless", "none"}:
        raise ValueError(
            f"{alpha_info.path} must be dimensionless; found BUNIT={alpha_info.unit!r}"
        )
    return {
        "stokes_i": str(stokes_info.path),
        "alpha": str(alpha_info.path),
        "spatial_shape": list(stokes_info.spatial_shape),
        "unit": stokes_info.unit,
        "reference_frequency_hz": entry.reference_frequency_hz,
    }


def write_image_model_previews(
    ctx: RunContext,
    center: SkyCoord,
    fov: u.Quantity,
) -> None:
    """Write FITS model previews for the weblog sky-model section."""
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    recenter = (center.ra.deg, center.dec.deg, fov.to(u.deg).value)
    for index, entry in enumerate(entries, start=1):
        image_path = primary_model_fits_path(entry)
        if image_path is None:
            continue
        suffix = "" if len(entries) == 1 else f"_{index:02d}"
        png_name = f"{ctx.work_dir.name}_fits_model{suffix}.png"
        png_path = ctx.work_dir / png_name
        write_fits_preview(
            image_path,
            png_path,
            "FITS Model",
            recenter=recenter,
            scale_factor=1000.0,
            bunit="mJy/pixel",
            colorbar_label="mJy/pixel",
        )
        ctx.manifest.add_output(
            "plot",
            png_name,
            role="fits_model",
            metadata={
                "model_entry_index": index - 1,
                "model_type": entry.type,
                "source_fits": str(image_path),
            },
        )


def inject_image_models(ctx: RunContext, visibility_path: Path) -> None:
    """Inject configured image models into an existing Measurement Set."""
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    ctx.add_milestone(
        "image_injection_started",
        "started",
        details={"n_model_entries": len(entries), "backend": "casa_ft"},
    )

    for index, entry in enumerate(entries):
        if isinstance(entry, StaticStokesMapsModelEntry):
            raise NotImplementedError(
                "static_stokes_maps is schema-ready, but the CASA backend path is "
                "planned for the next implementation phase."
            )
        if not isinstance(entry, ContinuumIAlphaModelEntry):
            continue
        report = validate_continuum_i_alpha(entry)
        product = prepare_continuum_i_alpha_for_casa(ctx, entry, index)
        run_casa_ft(
            visibility_path=visibility_path,
            model_paths=product.model_paths,
            nterms=product.nterms,
            reffreq=product.reffreq,
            incremental=index > 0,
        )
        ctx.manifest.add_output(
            "sky_model",
            product.model_paths[0].name,
            role="casa_model_image",
            metadata={
                "model_entry_index": index,
                "model_type": entry.type,
                "nterms": product.nterms,
                "reffreq": product.reffreq,
                "all_model_paths": [path.name for path in product.model_paths],
            },
        )
        ctx.add_milestone(
            "image_model_injected",
            "completed",
            details={
                "model_entry_index": index,
                "model_type": entry.type,
                "backend": "casa_ft",
                **report,
            },
        )

    merge_model_data_into_data(visibility_path)
    ctx.add_milestone(
        "image_injection_completed",
        "completed",
        details={"visibility_path": str(visibility_path), "model_data_merged": True},
    )


def prepare_continuum_i_alpha_for_casa(
    ctx: RunContext,
    entry: ContinuumIAlphaModelEntry,
    index: int,
) -> CasaModelProduct:
    """Create CASA image products for a continuum I+alpha model."""
    stokes_path = Path(entry.stokes_i).expanduser().resolve()
    alpha_path = Path(entry.alpha).expanduser().resolve()
    prefix = f"model_entry_{index + 1:02d}_continuum"
    tt0_fits = ctx.work_dir / f"{prefix}.tt0.fits"
    tt1_fits = ctx.work_dir / f"{prefix}.tt1.fits"
    tt0_image = ctx.work_dir / f"{prefix}.tt0.image"
    tt1_image = ctx.work_dir / f"{prefix}.tt1.image"

    shutil.copyfile(stokes_path, tt0_fits)
    with fits.open(stokes_path) as stokes_hdul, fits.open(alpha_path) as alpha_hdul:
        stokes_data = np.asarray(stokes_hdul[0].data, dtype=float)
        alpha_data = np.asarray(alpha_hdul[0].data, dtype=float)
        header = stokes_hdul[0].header.copy()
        header["BUNIT"] = stokes_hdul[0].header.get("BUNIT", "Jy/pixel")
        fits.writeto(tt1_fits, stokes_data * alpha_data, header=header, overwrite=True)

    for imagename in (tt0_image, tt1_image):
        if imagename.exists():
            shutil.rmtree(imagename)
    casa_tasks = import_casa_tasks()
    if casa_tasks is not None:
        importfits, _ = casa_tasks
        importfits(fitsimage=str(tt0_fits), imagename=str(tt0_image), overwrite=True)
        importfits(fitsimage=str(tt1_fits), imagename=str(tt1_image), overwrite=True)
    else:
        run_casa_importfits(
            ctx.work_dir,
            [(tt0_fits, tt0_image), (tt1_fits, tt1_image)],
        )

    return CasaModelProduct(
        model_paths=[tt0_image, tt1_image],
        nterms=2,
        reffreq=f"{entry.reference_frequency_hz}Hz",
        intermediates=[tt0_fits, tt1_fits],
    )


def import_casa_tasks():
    """Return in-process CASA tasks when they are importable in this Python env."""
    try:
        from casatasks import ft, importfits
    except Exception:
        return None
    return importfits, ft


def require_casa_tasks():
    """Import CASA tasks lazily and provide a clear runtime error if unavailable."""
    casa_tasks = import_casa_tasks()
    if casa_tasks is None:
        raise RuntimeError(
            "CASA casatasks.importfits and casatasks.ft are required for "
            "in-process image-model injection. Install casatasks in this "
            "environment or make the CASA executable available on PATH for "
            "batch-mode fallback."
        )
    return casa_tasks


def require_casa_executable() -> Path:
    """Return a CASA executable for batch-mode fallback."""
    executable = shutil.which("casa")
    if executable is None:
        raise RuntimeError(
            "CASA image-model injection requires either importable casatasks "
            "or a casa executable on PATH."
        )
    return Path(executable)


def run_casa_importfits(
    work_dir: Path,
    images: list[tuple[Path, Path]],
) -> None:
    """Run CASA importfits in batch mode for prepared FITS images."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_importfits.py"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import importfits",
        "except Exception:",
        "    pass",
    ]
    for fitsimage, imagename in images:
        lines.append(
            "importfits(fitsimage={!r}, imagename={!r}, overwrite=True)".format(
                str(fitsimage),
                str(imagename),
            )
        )
    run_casa_script(executable, script_path, lines)


def run_casa_ft(
    visibility_path: Path,
    model_paths: list[Path],
    nterms: int,
    reffreq: str,
    incremental: bool,
) -> None:
    """Run CASA ft into MODEL_DATA for one prepared model entry."""
    logger.info(
        f"CASA ft model={[str(path) for path in model_paths]} "
        f"nterms={nterms} reffreq={reffreq} incremental={incremental}"
    )
    casa_tasks = import_casa_tasks()
    if casa_tasks is not None:
        _, ft = casa_tasks
        ft(
            vis=str(visibility_path),
            model=[str(path) for path in model_paths],
            nterms=nterms,
            reffreq=reffreq,
            incremental=incremental,
            usescratch=True,
        )
        return

    executable = require_casa_executable()
    script_path = visibility_path.parent / "skasim_casa_ft.py"
    model_literal = "[" + ", ".join(repr(str(path)) for path in model_paths) + "]"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import ft",
        "except Exception:",
        "    pass",
        "ft(",
        f"    vis={str(visibility_path)!r},",
        f"    model={model_literal},",
        f"    nterms={int(nterms)},",
        f"    reffreq={reffreq!r},",
        f"    incremental={bool(incremental)!r},",
        "    usescratch=True,",
        ")",
    ]
    run_casa_script(executable, script_path, lines)


def run_casa_script(executable: Path, script_path: Path, lines: list[str]) -> None:
    """Write and execute one CASA batch script, surfacing useful failure output."""
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    command = [
        str(executable),
        "--nologger",
        "--nogui",
        "--log2term",
        "-c",
        str(script_path),
    ]
    logger.info(f"CASA batch command: {' '.join(command)}")
    result = subprocess.run(
        command,
        cwd=str(script_path.parent),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-80:])
        raise RuntimeError(
            f"CASA batch command failed with exit code {result.returncode}: "
            f"{script_path}\n{tail}"
        )


def merge_model_data_into_data(visibility_path: Path) -> None:
    """Add image-model MODEL_DATA into the delivered DATA column."""
    try:
        from casacore.tables import table
    except Exception as exc:
        raise RuntimeError(
            "python-casacore is required to merge MODEL_DATA into DATA."
        ) from exc

    with table(str(visibility_path), readonly=False, ack=False) as ms_table:
        columns = set(ms_table.colnames())
        if "DATA" not in columns or "MODEL_DATA" not in columns:
            raise ValueError(
                f"{visibility_path} must contain DATA and MODEL_DATA columns "
                "after CASA ft injection."
            )
        data = ms_table.getcol("DATA")
        model_data = ms_table.getcol("MODEL_DATA")
        ms_table.putcol("DATA", data + model_data)
