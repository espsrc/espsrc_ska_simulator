# Image-Model Smoke Configs

These configs exercise `continuum_i_alpha` image-model injection with compact
but useful observation settings: 60 seconds and 128 channels. They use the VLA
model FITS files in `models_for_testing/` and write outputs under
`demo_output/`.

The quick FITS configs use
`models_for_testing/CY4223_L_004_20180322_avg_target_3C320_vla_loop_01-0003-model.fits`.
That is a numbered WSClean per-channel model FITS, not the
`CY4223_L_004_20180322_avg_target_3C320_vla_loop_01-MFS-model.fits` product.
The generated flat-alpha map turns that one Stokes-I model plane into a
two-term CASA Taylor model for injection.

## Run all demos

From the repository root:

```bash
bash configs/fits_model/run_fits_model_demos.sh
```

## Prepare the alpha map

Run this once from the repository root:

```bash
PYTHONPATH=src python scripts/prepare_fits_model_alpha.py \
  models_for_testing/CY4223_L_004_20180322_avg_target_3C320_vla_loop_01-0003-model.fits \
  configs/fits_model/generated/3c320_flat_alpha.fits \
  --value -0.7 \
  --overwrite
```

## Quick full-execution checks

Use the environment that has Karabo, CASA, OSKAR, and optional WSClean installed.

```bash
conda run -n skasim skasim --config configs/fits_model/vla_c_image_only_quick.json
conda run -n skasim skasim --config configs/fits_model/meerkat_image_only_quick.json
conda run -n skasim skasim --config configs/fits_model/ska1mid_image_only_quick.json
```

Composite small JSON catalog plus FITS image model:

```bash
conda run -n skasim skasim --config configs/fits_model/vla_c_catalog_plus_image_quick.json
```

Optional WSClean image-product check:

```bash
conda run -n skasim skasim --config configs/fits_model/vla_c_image_only_wsclean.json
```

Inspect each run at:

```bash
demo_output/<run-name>/weblog.html
```

The Sky Model section should show a FITS Model preview, and the Science Products
section should show the downstream image product.

## Supported Image-Model Modes

### `continuum_i_alpha`

Use this for FITS images when you have:

```json
{
  "type": "continuum_i_alpha",
  "stokes_i": "path/to/stokes_i.fits",
  "alpha": "path/to/alpha.fits",
  "reference_frequency_hz": 1455500000.0
}
```

`stokes_i` must be a Jy/pixel-compatible FITS image. `alpha` must be a
dimensionless FITS image on the same spatial grid and WCS. `skasim` creates
CASA `.tt0` and `.tt1` images in the run directory and injects them with CASA
`ft`.

### `casa_taylor_terms`

Use this for existing CASA image-table Taylor terms, such as WSClean/CASA
products named `*.model.tt0` and `*.model.tt1`:

```json
{
  "type": "casa_taylor_terms",
  "tt0": "models_for_testing/Target_LSR_J1835_EVLA_L_inf_1_post.model.tt0",
  "tt1": "models_for_testing/Target_LSR_J1835_EVLA_L_inf_1_post.model.tt1",
  "reference_frequency_hz": 1575010361.7577474
}
```

This mode does not convert the pixel values through FITS for injection. `skasim`
copies the CASA image directories into the run directory, aligns their
single-channel spectral coordinate to `reference_frequency_hz`, and passes the
run-local copies to CASA `ft`. The `tt0` image is also exported to FITS for the
weblog model preview.

### `static_stokes_maps`

This schema accepts one or more FITS Stokes maps:

```json
{
  "type": "static_stokes_maps",
  "stokes_i": "path/to/i.fits",
  "stokes_q": "path/to/q.fits",
  "stokes_u": "path/to/u.fits",
  "stokes_v": "path/to/v.fits"
}
```

The schema is present, but backend injection is still reserved for the next
implementation phase.

## Long CASA Taylor-Term Examples

The configs under `configs/fits_model/casa_taylor_terms/` are examples for
large local model fixtures that are not generally distributed with the project:

- `Target_LSR_J1835_EVLA_L_inf_1_post.model.tt0/.tt1`
- `Target_M31STAR_EVLA_C_inf_1.model.tt0/.tt1`

Each pair has configs for:

- LOFAR HBA-like setup: 24 MHz centered at 144 MHz.
- VLA A L band: 1-2 GHz, represented as 1500 MHz center and 1000 MHz bandwidth.
- MeerKAT L band: 856-1712 MHz, represented as 1284 MHz center and 856 MHz bandwidth.
- SKA1-Mid Band 1b example: 700-1050 MHz, represented as 875 MHz center and 350 MHz bandwidth.

All of these long examples use 1024 channels and 600 seconds:

```bash
bash configs/fits_model/run_casa_taylor_model_demos.sh
```

These are intentionally heavier than the quick smoke configs. Their outputs are
written under `demo_output/casa_taylor_*`.

If an older run shows a blank LOFAR dirty image, check `skasim_casa_ft.log` for
`No overlap in frequency between image channels and selected data`. That means
the run used the original CASA image spectral coordinate rather than the
run-local aligned copy and should be rerun with the current code.
