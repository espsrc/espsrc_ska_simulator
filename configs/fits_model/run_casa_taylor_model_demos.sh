#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

for config in \
  configs/fits_model/casa_taylor_terms/lsr_j1835_lofar_144mhz.json \
  configs/fits_model/casa_taylor_terms/lsr_j1835_vla_a_lband.json \
  configs/fits_model/casa_taylor_terms/lsr_j1835_meerkat_lband.json \
  configs/fits_model/casa_taylor_terms/lsr_j1835_ska1mid_band1b.json \
  configs/fits_model/casa_taylor_terms/m31star_lofar_144mhz.json \
  configs/fits_model/casa_taylor_terms/m31star_vla_a_lband.json \
  configs/fits_model/casa_taylor_terms/m31star_meerkat_lband.json \
  configs/fits_model/casa_taylor_terms/m31star_ska1mid_band1b.json
do
  echo "=== Running ${config} ==="
  skasim --config "${config}"
done

echo "=== CASA Taylor-term model demo runs complete ==="
