# FITS Image-Model Smoke Configs

These configs exercise `continuum_i_alpha` image-model injection with the
smallest practical observation settings. They use the VLA model FITS files in
`models_for_testing/` and write outputs under `demo_output/`.

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

Composite catalog plus FITS image model:

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
