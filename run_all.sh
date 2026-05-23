#!/bin/bash
# run_all.sh — Run SKA Simulator pipelines using the custom 100 Gaussian source catalog.

# Exit immediately if any command fails
set -e

echo "=== SKASim Automated Test Suite ==="

# Create demo output directory
mkdir -p demo_output

# 0. Generate the deterministic reference catalog inside demo_output
# ==============================================================================
# SECTION: Configuration variables for all test runs
# ==============================================================================
SEED=42
FOV="1deg"
PIXELS=2048
CLEAN_ITERATIONS=10000
OBSERVATION_TIME=600
BANDWIDTH_MHZ=500
N_CHANNELS=2000
AUTO_THRESHOLD="1"
WSCLEAN_CMD="wsclean -weight briggs -0.5 -multiscale -auto-threshold ${AUTO_THRESHOLD} -channels-out 4"

echo "Generating a multi source reference Gaussian catalog with seed ${SEED}..."
conda run -n skasim python scripts/generate_gaussian_catalog.py --seed "${SEED}"

# ==============================================================================
# SECTION: Reference Gaussian Catalog WSClean Runs (demo_output)
# ==============================================================================
echo -e "\n--- Running Reference Gaussian Catalog with WSClean imager ---"

echo "Running Reference Catalog with MeerKAT..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope MeerKAT \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test1_reference_meerkat \
  --overwrite

echo "Running Reference Catalog with SKA-LOW-AAstar..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope SKA-LOW-AAstar \
  --telescope-version SKA_OST_ARRAY_CONFIG_2_3_1 \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test2_reference_ska_low_aastar \
  --overwrite

echo "Running Reference Catalog with SKA-MID-AA4..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope SKA-MID-AA4 \
  --telescope-version SKA_OST_ARRAY_CONFIG_2_3_1 \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test3_reference_ska_mid_aa4 \
  --overwrite

echo "Running Reference Catalog with SKA-MID-AAstar..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope SKA-MID-AAstar \
  --telescope-version SKA_OST_ARRAY_CONFIG_2_3_1 \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test4_reference_ska_mid_aastar \
  --overwrite

echo "Running Reference Catalog with SKA1LOW..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope SKA1LOW \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test5_reference_ska1low \
  --overwrite

echo "Running Reference Catalog with LOFAR..."
conda run -n skasim skasim \
  --model demo_output/reference_gaussian_catalog.json \
  --telescope LOFAR \
  --imager wsclean \
  --fov "${FOV}" \
  --pixels "${PIXELS}" \
  --clean-iterations "${CLEAN_ITERATIONS}" \
  --observation-time "${OBSERVATION_TIME}" \
  --bandwidth-mhz "${BANDWIDTH_MHZ}" \
  --n-channels "${N_CHANNELS}" \
  --wsclean-command "${WSCLEAN_CMD}" \
  --output-dir demo_output/test6_reference_lofar \
  --overwrite

echo -e "\n=== All tests defined in run_all.sh finished successfully! ==="



## Test 1: Generate a single random source
#conda run -n skasim skasim \
#  --output-dir demo_output/test1_single_source \
#  --catalog MIGHTEE \
#  --telescope MeerKAT \
#  --observation-time 30 \
#  --frequency-mhz 1300 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 2: Generate multiple random sources via command line
#conda run -n skasim skasim \
#  --output-dir demo_output/test2_multiple_sources \
#  --catalog MIGHTEE \
#  --telescope SKA-MID-AA4 \
#  --observation-time 30 \
#  --frequency-mhz 2000 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 3: Multiple sources with full Stokes parameters
#conda run -n skasim skasim \
#  --output-dir demo_output/test3_stokes_sources \
#  --catalog MIGHTEE \
#  --telescope LOFAR \
#  --observation-time 30 \
#  --frequency-mhz 144 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 4: Use built-in catalog (MIGHTEE)
#conda run -n skasim skasim \
#  --output-dir demo_output/test4_catalog_mightee \
#  --catalog MIGHTEE \
#  --telescope ASKAP \
#  --observation-time 30 \
#  --frequency-mhz 8000 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 5: Different telescope - SKA1MID
#conda run -n skasim skasim \
#  --output-dir demo_output/test5_ska1mid \
#  --catalog MIGHTEE \
#  --telescope SKA1MID \
#  --observation-time 30 \
#  --frequency-mhz 10000 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 6: Different telescope - VLA
#conda run -n skasim skasim \
#  --output-dir demo_output/test6_vla \
#  --catalog MIGHTEE \
#  --telescope VLA \
#  --observation-time 30 \
#  --frequency-mhz 5000 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 7: Different telescope - SKA1LOW (lower frequency)
#conda run -n skasim skasim \
#  --output-dir demo_output/test7_ska1low \
#  --catalog MIGHTEE \
#  --telescope SKA1LOW \
#  --observation-time 30 \
#  --frequency-mhz 150 \
#  --bandwidth-mhz 50 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
#
## Test 8: Test with WSClean imager instead of oskar-dirty
#conda run -n skasim skasim \
#  --output-dir demo_output/test8_wsclean \
#  --catalog MIGHTEE \
#  --telescope SKA-MID-AAstar \
#  --observation-time 30 \
#  --frequency-mhz 1300 \
#  --bandwidth-mhz 25 \
#  --n-channels 8 \
#  --pixels 256 \
#  --imager wsclean \
#  --clean-iterations 200 \
#  --overwrite
