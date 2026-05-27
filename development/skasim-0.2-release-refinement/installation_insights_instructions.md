# SKA Simulator 0.2 — Conda Installation Insights & Instructions

During the skasim 0.2 release refinement, the collaborator handoff noted that full simulator runs were blocked due to conda dependency solving failures for `karabo-pipeline`. Specifically, the environment solver failed to resolve transitive dependencies such as `fftw * mpi_mpich*` and `h5py * mpi_mpich*`.

This document captures the root cause analysis, outlines a robust and direct installation workflow, and provides verification instructions for future development.

---

## 1. Root Cause Analysis

### Symptom A: Conda solver failure (MPI packages)
When installing `karabo-pipeline` from the environment specification using default conda configuration, the solver threw unsatisfiable package conflict errors:
```text
LibMambaUnsatisfiableError: Encountered problems while solving:
  - package karabo-pipeline-0.33.2-py39hfeaa757_0 requires fftw * mpi_mpich*, but none of the providers can be installed
```

### Symptom B: OSKAR library loading failure (CUDA version mismatch)
If installing with `--no-channel-priority` without pinning `cuda-version`, the simulation failed at runtime with:
```text
RuntimeError: OSKAR library not found.
```
Inspecting the compiled `liboskar.so` with `ldd` revealed missing dynamic library links:
```text
libcudart.so.11.0 => not found
libcufft.so.10 => not found
```
This was caused by the environment installing newer `cuda-cudart` (e.g. `13.2.75`) and `libcufft` from `conda-forge` because of `--no-channel-priority`, whereas `oskar` in `i4ds` was compiled against CUDA 11 (`libcudart.so.11.0` and `libcufft.so.10`).

---

### Diagnosis & Solutions
1. **Disabling Strict Channel Priority (Fixes Symptom A)**: Bypassing strict channel priority with the `--no-channel-priority` flag allows the solver to search globally across all specified channels to fetch the compatible `mpi_mpich` builds of `fftw` and `h5py` from `conda-forge`.
2. **Explicitly Pinning the CUDA Version (Fixes Symptom B)**: When disabling strict channel priority, we must explicitly pin `cuda-version=11.7` in the install command. This prevents conda-forge from upgrading CUDA libraries to version 13, ensuring that the older, binary-compatible CUDA 11 runtime files required by `oskar` are preserved.


---

## 2. Step-by-Step Installation Instructions

Follow this exact procedure to initialize a fully functional `skasim` simulation runtime environment on this system. All steps run locally and install all dependencies fully self-contained in the conda environment.

### Step 1: Create the Clean Base Environment
Initialize the environment named `skasim` pinned to Python 3.9:
```bash
conda create -y -n skasim python=3.9
```

### Step 2: Install Karabo-Pipeline and Pin CUDA Runtime
Install the full radio astronomy runtime stack into the environment. We must relax the channel priority to find transitive dependencies from `conda-forge`, but explicitly constrain the CUDA runtime to `11.7` to prevent the solver from pulling incompatible newer CUDA versions from conda-forge:
```bash
conda install -y -n skasim --no-channel-priority \
  -c nvidia/label/cuda-11.7.0 \
  -c i4ds \
  -c conda-forge \
  karabo-pipeline "cuda-version=11.7"
```
*Note: This command cleanly solves the dependency tree, fetching the optimized `0.34.0` build of `karabo-pipeline` under Python 3.9 compiled against `libcudart.so.11.0` and `libcufft.so.10` alongside `mpich`, `casacore`, `everybeam`, and `wsclean` packages.*


### Step 3: Install the Local Repository in Editable Mode
Install the local `skasim` simulator package inside the environment in editable development mode:
```bash
conda run -n skasim pip install -e .
```

---

## 3. Environment Verification

Confirm that the runtime stack is functioning as expected by running the following tests:

### A. Print Command-Line Help
Check that the entry points are correctly bound and run without any pre-import issues:
```bash
conda run -n skasim skasim --help
```
*Expected output: Full `skasim` CLI arguments and parameters listing.*

### B. Verify OSKAR and WSClean Runtime Binaries
Check that the OSKAR Python module, OSKAR simulator binary, and WSClean binary
are present in the environment:
```bash
conda run -n skasim python -c "import karabo, oskar; print('karabo', karabo.__version__); print('oskar ok')"
conda run -n skasim which oskar_sim_interferometer
conda run -n skasim which wsclean
conda run -n skasim env OPENBLAS_NUM_THREADS=1 wsclean --version
```
*Expected output: `import oskar` succeeds, `oskar_sim_interferometer` resolves
inside the `skasim` environment, and `wsclean` resolves inside the same
environment. Direct WSClean commands need `OPENBLAS_NUM_THREADS=1`; the
`skasim` pipeline sets that automatically before launching WSClean, which is
why pipeline runs can work even if a bare `wsclean --version` aborts. `oskarpy`
is not required for this execution path.*

### C. List Simulated Telescopes
Test the integration with the `ska-sdp` and `karabo` backend libraries by listing simulated telescope configurations:
```bash
conda run -n skasim skasim --show-telescopes
```
*Expected output: A plain-text list of telescopes including `ALMA`, `VLA`, `MeerKAT`, `LOFAR`, and various `SKA` configurations.*

### D. Run a Complete WSClean Smoke Simulation
This command was verified on 2026-05-23. It exercises Karabo, OSKAR visibility
simulation, WSClean imaging, manifest writing, PNG preview generation, and the
weblog:
```bash
conda run -n skasim skasim \
  --output-dir smoke_mightee_wsclean_quick \
  --telescope MeerKAT \
  --observation-time 30 \
  --frequency-mhz 1300 \
  --bandwidth-mhz 25 \
  --n-channels 2 \
  --pixels 256 \
  --catalog MIGHTEE \
  --imager wsclean \
  --clean-iterations 20 \
  --overwrite
```
*Expected output: exit code 0, manifest `status: completed`, a
`visibilities.MS` directory, WSClean FITS/PNG image products, and
`weblog.html` in `smoke_mightee_wsclean_quick/`.*

### E. Run the Pytest Test Suite
Validate the entire simulator codebase against the complete test suite:
```bash
conda run -n skasim pytest
```
*Expected output: the normal unit-test suite passes. The lightweight missing-
Karabo test is intended for a pip-only environment and should be interpreted
accordingly when run inside the full `skasim` runtime environment.*

---

## 4. Key Takeaways for Future Handoffs

- **Solver Strategy**: In complex scientific stacks with nested C/C++ libraries (e.g. `casacore`, `everybeam`, `wsclean`), channel priority conflicts are very common. Always consider installing with `--no-channel-priority` if conflicts arise.
- **Python Compatibility**: Pinned dependency requirements restrict version increments. For the `skasim` 0.2 series, **Python 3.9** remains the gold standard, ensuring full compatibility with both older and newer `karabo-pipeline` channels.
- **OSKAR Runtime**: OSKAR is installed in the working environment. The
  supported verification is `import oskar`, `which oskar_sim_interferometer`,
  and a completed visibility simulation. `oskarpy` is not part of the required
  runtime path.
