"""Public façade for image-model loading, validation and injection."""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from loguru import logger

from ...config import (
    CasaTaylorTermsModelEntry,
    ComponentSkyModelEntry,
    ContinuumIAlphaModelEntry,
    ImgConfig,
    ModelEntry,
    SimConfig,
    SpectralCubeModelEntry,
    StaticStokesMapsModelEntry,
    has_spectral_cube_model,
    spectral_cube_model_entries,
)
from ...manifest import RunContext
import subprocess

from ...runtime import require_casacore as require_casacore
from .casa_interop import (
    CasaModelProduct,
    _resample_spectral_axis_to_ms_channels,
    _set_crval4_via_script,
    adjust_spectral_reference,
    merge_model_data_into_data,
    prepare_casa_taylor_terms,
    prepare_continuum_i_alpha_for_casa,
    prepare_spectral_cube_for_casa,
    require_casa_executable,
    run_casa_exportfits,
    run_casa_ft,
    run_casa_importfits,
    run_casa_set_spectral_coordinate,
    validate_casa_taylor_terms,
)
from .fits_io import (
    FitsCubeInfo,
    FitsImageInfo,
    _select_wsclean_img_config,
    component_model_entries,
    image_model_center,
    image_model_entries,
    primary_model_fits_path,
    read_fits_cube_info,
    read_fits_image_info,
    validate_continuum_i_alpha,
    validate_spectral_cube,
)
from .previews import (
    run_moment8_for_spectral_cube,
    write_image_model_previews,
    write_spectral_cube_input_preview,
)


def inject_image_models(ctx: RunContext, visibility_path: Path) -> None:
    """Inject configured image models into an existing Measurement Set."""
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    backends = {"casa_ft", "wsclean_predict"}
    ctx.add_milestone(
        "image_injection_started",
        "started",
        details={"n_model_entries": len(entries), "backends": sorted(backends)},
    )

    for index, entry in enumerate(entries):
        if isinstance(entry, StaticStokesMapsModelEntry):
            raise NotImplementedError(
                "static_stokes_maps is schema-ready, but the CASA backend path is "
                "planned for the next implementation phase."
            )
        if isinstance(entry, ContinuumIAlphaModelEntry):
            img_config = _select_wsclean_img_config(ctx.config.imaging)
            report = validate_continuum_i_alpha(entry)
            product = prepare_continuum_i_alpha_for_casa(ctx, entry, index)

            if entry.injection_backend == "casa_ft":
                logger.warning(
                    "continuum_i_alpha injection_backend='casa_ft' is deprecated; "
                    "prefer 'wsclean_predict'."
                )
                run_casa_ft(
                    visibility_path=visibility_path,
                    model_paths=product.model_paths,
                    nterms=product.nterms,
                    reffreq=product.reffreq,
                    incremental=index > 0,
                )
                backend = "casa_ft"
            else:
                from ..wsclean_predict import inject_continuum_i_alpha_with_wsclean_predict

                report_predict = inject_continuum_i_alpha_with_wsclean_predict(
                    ctx,
                    entry,
                    index,
                    visibility_path,
                    img_config,
                    product,
                )
                backend = "wsclean_predict"
                report = {**report, **report_predict}
        elif isinstance(entry, CasaTaylorTermsModelEntry):
            report = validate_casa_taylor_terms(entry)
            product = prepare_casa_taylor_terms(ctx, entry, index)
            logger.warning(
                "casa_taylor_terms uses the deprecated CASA ft backend; "
                "consider migrating to continuum_i_alpha with wsclean_predict."
            )
            run_casa_ft(
                visibility_path=visibility_path,
                model_paths=product.model_paths,
                nterms=product.nterms,
                reffreq=product.reffreq,
                incremental=index > 0,
            )
            backend = "casa_ft"
        elif isinstance(entry, SpectralCubeModelEntry):
            img_config = _select_wsclean_img_config(ctx.config.imaging)
            report = validate_spectral_cube(entry, ctx.config.observation, img_config)
            product = prepare_spectral_cube_for_casa(ctx, entry, index, report)
            from ..wsclean_predict import inject_spectral_cube_with_wsclean_predict

            assert product.cube_data is not None
            assert product.header is not None
            assert product.freq_axis is not None
            report_predict = inject_spectral_cube_with_wsclean_predict(
                ctx,
                entry,
                index,
                visibility_path,
                img_config,
                product.cube_data,
                product.header,
                product.freq_axis,
            )
            backend = "wsclean_predict"
            report = {**report, **report_predict}
        else:
            continue
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
                "backend": backend,
                **report,
            },
        )

    merge_model_data_into_data(visibility_path)
    ctx.add_milestone(
        "image_injection_completed",
        "completed",
        details={"visibility_path": str(visibility_path), "model_data_merged": True},
    )


__all__ = [
    "CasaModelProduct",
    "FitsCubeInfo",
    "FitsImageInfo",
    "adjust_spectral_reference",
    "component_model_entries",
    "has_spectral_cube_model",
    "image_model_center",
    "image_model_entries",
    "inject_image_models",
    "merge_model_data_into_data",
    "prepare_casa_taylor_terms",
    "prepare_continuum_i_alpha_for_casa",
    "prepare_spectral_cube_for_casa",
    "primary_model_fits_path",
    "read_fits_cube_info",
    "read_fits_image_info",
    "require_casa_executable",
    "run_casa_exportfits",
    "run_casa_ft",
    "run_casa_importfits",
    "run_casa_set_spectral_coordinate",
    "run_moment8_for_spectral_cube",
    "spectral_cube_model_entries",
    "validate_casa_taylor_terms",
    "validate_continuum_i_alpha",
    "validate_spectral_cube",
    "write_image_model_previews",
    "write_spectral_cube_input_preview",
]
