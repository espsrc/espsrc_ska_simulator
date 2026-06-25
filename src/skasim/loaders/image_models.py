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

from ..config import (
    CasaTaylorTermsModelEntry,
    ComponentSkyModelEntry,
    ContinuumIAlphaModelEntry,
    ModelEntry,
    SimConfig,
    StaticStokesMapsModelEntry,
)
from ..imaging import write_fits_preview
from ..manifest import RunContext
from ..runtime import require_casacore


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
        entry for entry in config.models if isinstance(entry, ComponentSkyModelEntry)
    ]


def image_model_entries(config: SimConfig) -> list[ModelEntry]:
    return [
        entry
        for entry in config.models
        if isinstance(
            entry,
            (
                ContinuumIAlphaModelEntry,
                CasaTaylorTermsModelEntry,
                StaticStokesMapsModelEntry,
            ),
        )
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

    for index, entry in enumerate(entries, start=1):
        image_path = primary_model_fits_path(entry)
        export_path = None
        if image_path is None and isinstance(entry, CasaTaylorTermsModelEntry):
            image_path = Path(entry.tt0).expanduser().resolve()
            export_path = ctx.work_dir / f"model_entry_{index:02d}_casa_taylor.tt0.fits"
        if image_path is None:
            continue
        # use the FITS image's own WCS center to avoid recentering NaN
        # when the model and sky-catalog coordinates differ
        try:
            info = read_fits_image_info(image_path)
            if info.center is None:
                recenter = None
            else:
                assert isinstance(info.center, SkyCoord)  # narrow for type checker
                recenter = (
                    info.center.ra.deg,
                    info.center.dec.deg,
                    fov.to(u.deg).value,
                )
        except Exception:
            recenter = None
        suffix = "" if len(entries) == 1 else f"_{index:02d}"
        png_name = f"{ctx.work_dir.name}_fits_model{suffix}.png"
        png_path = ctx.work_dir / png_name
        preview_source = image_path
        if export_path is not None:
            run_casa_exportfits(ctx.work_dir, image_path, export_path)
            preview_source = export_path
        write_fits_preview(
            preview_source,
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
                "preview_fits": str(preview_source),
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
        if isinstance(entry, ContinuumIAlphaModelEntry):
            report = validate_continuum_i_alpha(entry)
            product = prepare_continuum_i_alpha_for_casa(ctx, entry, index)
        elif isinstance(entry, CasaTaylorTermsModelEntry):
            report = validate_casa_taylor_terms(entry)
            product = prepare_casa_taylor_terms(ctx, entry, index)
        else:
            continue
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


def validate_casa_taylor_terms(entry: CasaTaylorTermsModelEntry) -> dict:
    """Validate an existing CASA Taylor-term image model entry."""
    model_paths = [
        Path(path).expanduser().resolve()
        for path in (entry.tt0, entry.tt1)
        if path is not None
    ]
    if not model_paths:
        raise ValueError("casa_taylor_terms requires at least tt0.")
    for path in model_paths:
        if not path.is_dir():
            raise ValueError(f"{path} must be a CASA image table directory.")
        if not (path / "table.dat").exists():
            raise ValueError(f"{path} does not look like a CASA image table.")
    return {
        "model_paths": [str(path) for path in model_paths],
        "nterms": len(model_paths),
        "reference_frequency_hz": entry.reference_frequency_hz,
    }


def prepare_casa_taylor_terms(
    ctx: RunContext,
    entry: CasaTaylorTermsModelEntry,
    index: int,
) -> CasaModelProduct:
    """Copy CASA Taylor-term images into the run and align their spectral reference.

    The reference frequency is adjusted to the observation band centre.
    For nterms≥2, tt0 pixel data is scaled:  tt0' = tt0 · (ν_obs / ν_old)^α
    where α = mean(tt1) / mean(tt0).  tt1 pixel data is unchanged.
    For nterms=1, only CRVAL4 is updated.
    """
    new_ref_hz = ctx.config.observation.frequency_mhz * 1e6
    old_ref_hz = entry.reference_frequency_hz

    source_paths = [
        Path(path).expanduser().resolve()
        for path in (entry.tt0, entry.tt1)
        if path is not None
    ]
    prefix = f"model_entry_{index + 1:02d}_casa_taylor"
    model_paths = []
    for term_index, source_path in enumerate(source_paths):
        target_path = ctx.work_dir / f"{prefix}.tt{term_index}.image"
        if target_path.exists():
            shutil.rmtree(target_path)
        shutil.copytree(source_path, target_path)
        model_paths.append(target_path)

    nterms = len(model_paths)

    if nterms >= 2:
        # alpha map: element-wise tt1/tt0, with safeguard for tt0==0 (those pixels
        # stay zero regardless of spectral index).
        casacore_table = require_casacore()
        with casacore_table(str(model_paths[0]), readonly=True, ack=False) as tbl0:
            tt0_data = np.asarray(tbl0.getcol("map"))
        with casacore_table(str(model_paths[1]), readonly=True, ack=False) as tbl1:
            tt1_data = np.asarray(tbl1.getcol("map"))

        # compute element-wise alpha, safeguarding divide-by-zero
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha_map = np.where(
                tt0_data != 0,
                tt1_data / tt0_data,
                0.0,
            )
        alpha_mean = float(np.mean(alpha_map))
        logger.info(
            f"prepare_casa_taylor_terms: alpha_map_mean={alpha_mean:.6f} "
            f"from element-wise tt1/tt0"
        )

        adjust_spectral_reference(
            model_paths[0],
            old_ref_hz,
            new_ref_hz,
            alpha_map=alpha_map,
        )
        # tt1: only set CRVAL4, no pixel-data correction
        adjust_spectral_reference(
            model_paths[1],
            old_ref_hz,
            new_ref_hz,
            alpha_map=None,
        )
    else:
        # nterms=1: spectrally flat, only set CRVAL4
        adjust_spectral_reference(
            model_paths[0],
            old_ref_hz,
            new_ref_hz,
            alpha_map=None,
        )

    ctx.add_milestone(
        "adjusted_spectral_reference",
        "completed",
        details={
            "model_type": "casa_taylor_terms",
            "old_reference_frequency_hz": old_ref_hz,
            "new_reference_frequency_hz": new_ref_hz,
            "nterms": nterms,
        },
    )

    return CasaModelProduct(
        model_paths=model_paths,
        nterms=nterms,
        reffreq=f"{new_ref_hz}Hz",
        intermediates=model_paths,
    )


def prepare_continuum_i_alpha_for_casa(
    ctx: RunContext,
    entry: ContinuumIAlphaModelEntry,
    index: int,
) -> CasaModelProduct:
    """Create CASA image products for a continuum I+alpha model.

    Adjusts the spectral reference to the observation band centre using
    the explicit spectral index from the model entry.
    """
    new_ref_hz = ctx.config.observation.frequency_mhz * 1e6
    old_ref_hz = entry.reference_frequency_hz

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

    # adjust spectral reference: tt0 pixel data scaled by (ν_new/ν_old)^α
    # α is the mean spectral index from the alpha map
    alpha_mean = float(np.mean(alpha_data))
    adjust_spectral_reference(
        tt0_image,
        old_ref_hz,
        new_ref_hz,
        alpha_map=alpha_data,
    )
    # tt1: only set CRVAL4, no pixel-data correction
    adjust_spectral_reference(
        tt1_image,
        old_ref_hz,
        new_ref_hz,
        alpha_map=None,
    )

    ctx.add_milestone(
        "adjusted_spectral_reference",
        "completed",
        details={
            "model_type": "continuum_i_alpha",
            "old_reference_frequency_hz": old_ref_hz,
            "new_reference_frequency_hz": new_ref_hz,
            "alpha_mean": alpha_mean,
            "nterms": 2,
        },
    )

    return CasaModelProduct(
        model_paths=[tt0_image, tt1_image],
        nterms=2,
        reffreq=f"{new_ref_hz}Hz",
        intermediates=[tt0_fits, tt1_fits],
    )


def adjust_spectral_reference(
    image_path: Path,
    old_ref_hz: float,
    new_ref_hz: float,
    alpha_map: np.ndarray | None = None,
) -> float:
    """Adjust the spectral reference of a CASA image to the observation band centre.

    For nterms=1 (alpha_map is None), set CRVAL4 to new_ref_hz only — the model
    is spectrally flat and no pixel-data correction is needed.

    For nterms≥2 (alpha_map provided), correct the pixel data element-wise following
    CASA's Taylor-series convention::

        tt0'(x,y) = tt0(x,y) · (ν_new / ν_old) ^ α(x,y)

    where α(x,y) = tt1(x,y) / tt0(x,y)   for each pixel.  CRVAL4 is also set.

    Returns the adjusted reference frequency in Hz (always new_ref_hz).
    """
    if alpha_map is not None:
        casacore_table = require_casacore()

        # element-wise scaling
        factor = (new_ref_hz / old_ref_hz) ** alpha_map
        with casacore_table(str(image_path), readonly=False, ack=False) as tbl:
            data = np.asarray(tbl.getcol("map"))
            corrected = data * factor
            tbl.putcol("map", corrected)
        logger.info(
            f"Adjusted pixel data in {image_path}: ref_freq {old_ref_hz:.3e}Hz → "
            f"{new_ref_hz:.3e}Hz (factor range: [{np.min(factor):.4f}, {np.max(factor):.4f}])"
        )
    else:
        logger.info(
            f"Spectrally flat image {image_path}: updating CRVAL4 = {new_ref_hz:.3e}Hz "
            "(no pixel-data correction)"
        )

    _set_crval4_via_script(image_path.parent, [image_path], new_ref_hz)
    return new_ref_hz


def _image_has_spectral_axis(image_path: Path) -> bool:
    """Return True if the CASA image has a frequency/spectral axis."""
    casacore_table = require_casacore()
    with casacore_table(str(image_path), readonly=True, ack=False) as tbl:
        coords = tbl.getcolkeywords("map")
    dim_names = [str(name).lower() for name in coords.get("dimnames", [])]
    return "frequency" in dim_names


def _set_crval4_via_script(
    work_dir: Path,
    image_paths: list[Path],
    frequency_hz: float,
) -> None:
    """Set CRVAL4 on CASA images — subprocess fallback for environments without casatasks."""
    # 2D images have no spectral axis; CASA ft will treat them as spectrally flat
    # with the reference frequency passed as reffreq, so CRVAL4 cannot (and need not) be set.
    image_paths = [p for p in image_paths if _image_has_spectral_axis(p)]
    if not image_paths:
        logger.info(
            "Spectrally flat 2D image(s): skipping CRVAL4 update; ft reffreq carries the reference frequency"
        )
        return

    try:
        from casatasks import imhead

        for image_path in image_paths:
            imhead(
                imagename=str(image_path),
                mode="put",
                hdkey="crval4",
                hdvalue=f"{frequency_hz}Hz",
            )
        return
    except Exception:
        pass  # fall through to batch mode

    run_casa_set_spectral_coordinate(work_dir, image_paths, frequency_hz)


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
        "except Exception as _e:",
        "    raise RuntimeError('casatasks.importfits is not available: ' + str(_e))",
    ]
    for fitsimage, imagename in images:
        lines.append(
            "importfits(fitsimage={!r}, imagename={!r}, overwrite=True)".format(
                str(fitsimage),
                str(imagename),
            )
        )
    run_casa_script(executable, script_path, lines)
    # verify side-effects; if not present, casatasks ran but failed silently
    for _fitsimage, imagename in images:
        table_dat = imagename / "table.dat"
        if not table_dat.exists():
            raise RuntimeError(
                f"CASA importfits did not create {imagename} as expected."
            )


def run_casa_exportfits(
    work_dir: Path,
    imagename: Path,
    fitsimage: Path,
) -> None:
    """Run CASA exportfits in batch mode for a CASA image table."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_exportfits.py"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import exportfits",
        "except Exception:",
        "    pass",
        "exportfits(imagename={!r}, fitsimage={!r}, overwrite=True)".format(
            str(imagename),
            str(fitsimage),
        ),
    ]
    run_casa_script(executable, script_path, lines)


def run_casa_set_spectral_coordinate(
    work_dir: Path,
    image_paths: list[Path],
    frequency_hz: float,
) -> None:
    """Set the single-channel spectral coordinate of CASA images to the run reference."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_set_spectral_coordinate.py"
    hdvalue = f"{frequency_hz}Hz"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import imhead",
        "except Exception:",
        "    pass",
    ]
    for image_path in image_paths:
        lines.append(
            "imhead(imagename={!r}, mode='put', hdkey='crval4', hdvalue={!r})".format(
                str(image_path),
                hdvalue,
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
        "except Exception as _e:",
        "    raise RuntimeError('casatasks.ft is not available: ' + str(_e))",
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
    log_path = script_path.with_suffix(".log")
    command = [
        str(executable),
        "--nologger",
        "--nogui",
        "--log2term",
        "-c",
        str(script_path),
    ]
    logger.info(f"CASA batch command: {' '.join(command)}")
    with log_path.open("w", encoding="utf-8") as log_file:
        result = subprocess.run(
            command,
            cwd=str(script_path.parent),
            text=True,
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
        if result.returncode != 0:
            tail = "\n".join(log_path.read_text(encoding="utf-8").splitlines()[-80:])
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
