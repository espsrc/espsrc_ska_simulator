#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

ALPHA_PATH="configs/fits_model/generated/3c320_flat_alpha.fits"
STOKES_I="models_for_testing/CY4223_L_004_20180322_avg_target_3C320_vla_loop_01-0003-model.fits"

PYTHONPATH=src python scripts/prepare_fits_model_alpha.py \
  "${STOKES_I}" \
  "${ALPHA_PATH}" \
  --value -0.7 \
  --overwrite

for config in \
  configs/fits_model/vla_c_image_only_quick.json \
  configs/fits_model/meerkat_image_only_quick.json \
  configs/fits_model/ska1mid_image_only_quick.json \
  configs/fits_model/vla_c_catalog_plus_image_quick.json \
  configs/fits_model/vla_c_image_only_wsclean.json
do
  echo "=== Running ${config} ==="
  conda run -n skasim skasim --config "${config}"
done

echo "=== FITS model demo runs complete ==="
