# skasim 0.2 Release Refinement

This directory is the collaborator handoff for the skasim 0.2 release
refinement. Start here, then use the linked files for detail.

## Reading Order

- [PRD](./PRD.md): original motivation, scope, user stories, and acceptance
  criteria.
- [Issue index and dependency graph](./ISSUES.md): the 12 implementation slices.
- [Local issues](./issues/): per-slice acceptance criteria and suggested commit
  messages.
- [Assessment 01](./assessment_01.md): external review of the first 12 commits.
- [Assessment 01 resolution](./assessment_01_resolution.md): accepted review
  findings and final decisions.
- [Installation Insights Instructions](./installation_insights_instructions.md): step-by-step instructions to successfully install the Conda environment and troubleshoot Karabo dependency resolver issues.


## Current State

- Release target is `0.2.0`, reflected in `pyproject.toml`, `docs/conf.py`, and
  the public README badge.
- All 12 planned issues have been implemented and then refined after
  `assessment_01.md`.
- The current branch still has uncommitted post-assessment changes.
- Domain language should follow the repository `CONTEXT.md`: catalog, sky
  model, sky model source, source intensity, run, image product, WSClean
  command, manifest, and weblog.

## Why This Work Was Needed

The 0.1 interface mixed domain concepts and implementation details:

- Built-in catalogs were selected by numeric IDs, which made commands and
  manifests hard to interpret.
- `--I` was ambiguous: it looked like one Stokes component but actually created
  generated source intensities.
- `--cleaning` was a boolean switch for an imaging backend instead of an
  explicit image product choice.
- WSClean was hardcoded as `wsclean`, which did not fit systems where WSClean is
  provided through Singularity.
- Karabo imports happened too early for pip-only CLI/help/config workflows.
- Run records were not structured enough for reliable output inspection.
- Failed runs did not always produce an inspectable weblog.
- WSClean output discovery was fragile and could mix old outputs with current
  outputs.

The 0.2 work consolidates the interface around explicit domain language and
makes the runtime behavior easier to test and inspect.

## Implementation Timeline

Initial 12 commits:

- `bdeba01 docs: add skasim 0.2 development baseline`
- `cf2794b refactor: isolate Karabo behind runtime imports`
- `a4d4378 docs: add pip and conda installation paths`
- `6cc4920 feat: use named sky model sources`
- `365f2fe feat: clarify generated source intensity CLI`
- `3e32dac feat: replace cleaning flag with explicit imager`
- `a54f645 feat: configure WSClean command execution`
- `d7c22d3 refactor: make WSClean outputs run scoped`
- `183b30a feat: structure run manifest outputs`
- `3c74c55 feat: render weblogs for every run`
- `d71f689 refactor: deepen run records and run ids`
- `7a254ec docs: finalize skasim 0.2 release examples`

Post-assessment changes are currently uncommitted. They address
`assessment_01.md` and are summarized in
[assessment_01_resolution.md](./assessment_01_resolution.md).

## New Functionality

- Pip-only imports and `python -m skasim.cli --help` work without Karabo.
- Full simulation execution loads Karabo only at runtime and fails with a clear
  installation message when Karabo is unavailable.
- Built-in catalogs are named: `MIGHTEE`, `GLEAM`, and `SKAMid`.
- A run accepts one explicit sky model source: file-backed `--model`, named
  `--catalog`, or generated source flux densities.
- Generated source Stokes I flux densities use `--flux-density`.
- Optional generated source polarization uses `--stokes-q`, `--stokes-u`, and
  `--stokes-v`.
- Observing frequency, channel width, and duration use explicit unit-bearing
  flags such as `--frequency-mhz`, `--channel-width-mhz`, and
  `--observation-time`.
- `--output-dir` names the exact output directory.
- Imaging is selected explicitly with `--imager oskar-dirty` or
  `--imager wsclean`.
- `--list-telescopes` separates telescopes that need no version from those
  that require `--telescope-version`, and shows accepted version names.
- WSClean iterations use `--clean-iterations`.
- WSClean invocation is configurable through `--wsclean-command`.
- WSClean commands are parsed into argv and executed with `shell=False`.
- WSClean outputs use a run-scoped prefix and are recorded as structured image
  products.
- WSClean `-channels-out` is capped to the simulated channel count, up to 8,
  so reduced-channel smoke runs are valid.
- FITS image-product PNG previews are rendered with APLpy and CMasher using
  in-memory 2D celestial HDUs, a non-interactive Matplotlib backend, WCS axes,
  `cmr.rainforest`, mJy/beam color scales, RMS-based clipping, optional
  contours, and beam overlays when FITS beam metadata is available.
- The weblog groups WSClean MFS model, clean, and residual products as the
  primary science products. The PSF is linked from the product header rather
  than shown as a fourth image. Dirty images are shown only as a fallback when
  model/clean/residual products are not available, such as `oskar-dirty` runs.
- The weblog keeps telescope and observation metadata together, including the
  resolved frequency range, central frequency, channel width, total bandwidth,
  channel count, observation time, and timestep-derived integration time when
  available. Imaging setup is intentionally compact: image size, total FoV,
  pixel size, and representative FITS beam metadata when present. Cleaning and
  imager setup lists the effective WSClean or OSKAR parameters, including
  WSClean weighting, multiscale, mgain, auto-threshold, auto-mask, channels-out,
  join-channels, local-rms, output prefix, and visibility input. Science
  products display representative beam metadata once per product instead of
  repeating it under every preview image.
- The weblog ends with a compact software-version table for reproducibility:
  `skasim`, Karabo, OSKAR, WSClean, and Python. Python package metadata is used
  when available, and Conda metadata is used for runtime tools such as WSClean
  and OSKAR that are installed as Conda packages.
- Each run writes two sky-model preview plots after sky-model loading: the full
  source model and a FoV-matched zoom. The weblog shows these in the Sky Model
  section. Source shapes are rendered as ellipses from major axis, minor axis,
  and position angle metadata. The sky-model colormap is reversed so brighter
  sources render darker and faint sources render lighter. Source position
  angles are converted from astronomical PA, measured east of north, to the
  Matplotlib convention. Sources smaller than the plotted FoV-dependent
  resolution threshold are drawn as flux-scaled crosses so compact components
  remain visible in full-field plots.
- JSON sky-model loading preserves full source metadata, not only reduced
  `(ra, dec, I)` rows, so preview plots and downstream records retain source
  sizes, position angles, spectral metadata, and polarization fields.
- The reference Gaussian JSON catalog generator assigns per-source spectral
  indices from a normal distribution centered on `-0.5` with standard deviation
  `0.2`; JSON ingestion preserves these values in the resulting sky model. It
  now produces a broader demonstration population: most sources are centered
  around the 1 mJy scale, two deterministic bright sources are near 60-100 mJy,
  several faint sources are in the 10-20 microJy range, and source sizes span
  sub-arcsec compact components through several-arcmin extended Gaussians. The
  generator also writes a matching DS9 FK5 region file with one ellipse and one
  compact cross marker per source.
- The weblog places the telescope layout next to the compact
  observation/telescope configuration, keeps run timing small in the header,
  and moves the detailed milestone timeline after the science images.
- Every run records logs, manifest, visibility data, image products, plots, sky
  model source, and weblog outputs by structured kind where applicable.
- Every run writes `weblog.html`, including failed runs.
- Default run IDs use second precision: `YYYYMMDD_HHMMSS_<telescope>`.
- Existing output directories require `--overwrite`; there is no interactive
  prompt. With `--overwrite`, the exact output directory is removed and
  recreated before the run starts so stale MeasurementSet locks, WSClean
  products, logs, and manifests cannot mix with the new run.
- CLI and pipeline execution force a non-interactive Matplotlib backend and
  close telescope-layout figures after saving. This avoids X11/ICE shutdown
  errors on headless hosts after otherwise successful OSKAR/WSClean runs.

## Removed Or Rejected

- CLI `--I` is removed and fails with a migration message.
- CLI `--stokes-i`, `--fits`, `--json`, `--json-fg`, `--Q`,
  `--U`, `--V`, `--ref-freq`, `--freq`, `--bandwidth`, `--delta-freq`,
  `--seconds`, `--prefix`, `--niter`, `--scale-I`, and `--imaging-niter` are
  removed or renamed and fail with migration messages.
- CLI `--cleaning` is removed and fails with a migration message.
- Numeric catalog IDs such as `--catalog 1` are removed and fail with a
  migration message.
- Python config fields `SimConfig.I`, `SimConfig.Q`, `SimConfig.U`,
  `SimConfig.V`, `SimConfig.cleaning`, `SimConfig.source_names`,
  `SimConfig.ref_freq_hz`, `SimConfig.json_fg`, `SimConfig.output_prefix`,
  `SimConfig.niter`, and `SimConfig.scale_I` are removed.
- Python config field `ImgConfig.algorithm` is removed; `ImgConfig.imager` is
  the source of truth.
- Config models reject unknown fields with Pydantic `extra="forbid"`.
- Scientific/configuration metadata is no longer encoded in WSClean output
  filenames; it belongs in the manifest.
- The proposed full Python f-string weblog rewrite was rejected. The weblog
  remains Jinja2-based so HTML structure stays separated from Python data
  preparation.

## Modified Files By Area

- Release metadata and user docs: `pyproject.toml`, `README.md`,
  `docs/conf.py`, `docs/installation.rst`, `docs/guide.rst`,
  `docs/examples.rst`, `docs/introduction.rst`.
- Planning and handoff: `development/skasim-0.2-release-refinement/README.md`,
  `ISSUES.md`, `PRD.md`, `issues/*.md`, `assessment_01.md`,
  `assessment_01_resolution.md`.
- Public API and configuration: `src/skasim/__init__.py`,
  `src/skasim/config.py`.
- CLI and runtime loading: `src/skasim/cli.py`, `src/skasim/runtime.py`.
- Sky model loading: `src/skasim/sky.py`, `src/skasim/fits_helper.py`.
- Simulation orchestration: `src/skasim/pipeline.py`.
- Imaging: `src/skasim/imaging.py`.
- Run records and weblog: `src/skasim/manifest.py`, `src/skasim/weblog.py`,
  `src/skasim/templates/weblog.html.j2`.
- Tests: `tests/test_cli.py`, `tests/test_config.py`,
  `tests/test_imaging.py`, `tests/test_manifest.py`, `tests/test_pipeline.py`,
  `tests/test_release_docs.py`, `tests/test_runtime_imports.py`,
  `tests/test_weblog.py`.

## Runtime Notes

Lightweight verification uses a temporary dependency path:

```bash
PYTHONPATH=/tmp/skasim-pydeps:src python -m py_compile src/skasim/*.py
PYTHONPATH=/tmp/skasim-pydeps:src python -m pytest -q
```

Latest result:

```text
109 passed in the focused CLI/config/pipeline/manifest/imaging/runtime/docs/weblog subset
Reference MeerKAT JSON-catalog WSClean command completed with exit code 0 after clean overwrite/backend fixes.
```

The temporary dependency path is only for this development session. The verified
full-runtime installation path is:

```bash
conda create -y -n skasim python=3.9
conda install -y -n skasim --no-channel-priority \
  -c nvidia/label/cuda-11.7.0 \
  -c i4ds \
  -c conda-forge \
  karabo-pipeline "cuda-version=11.7"
conda run -n skasim pip install -e .
```

The Conda environment uses the supported `skasim` env name and Karabo conda path.
The initial handoff was blocked by Conda resolver issues with transitive
dependencies (e.g., `fftw * mpi_mpich*`) and then by an OSKAR CUDA runtime
mismatch. The working fix is both:

- disable strict channel priority so MPI dependencies can resolve;
- pin `cuda-version=11.7` so the OSKAR library finds `libcudart.so.11.0` and
  `libcufft.so.10`.

For exact commands to install and verify the runtime successfully, see
[installation_insights_instructions.md](./installation_insights_instructions.md).

Runtime verification commands:

```bash
conda run -n skasim python -c "import karabo, oskar; print('karabo', karabo.__version__); print('oskar ok')"
conda run -n skasim which oskar_sim_interferometer
conda run -n skasim which wsclean
conda run -n skasim skasim --help
```

`oskarpy` is not installed and is not required. The working path uses the
`oskar` module and the OSKAR backend through Karabo; the successful smoke runs
created `visibilities.MS`.


Canonical smoke run command shape after the CLI cleanup:

```bash
conda run -n skasim skasim \
  --output-dir smoke_mightee_wsclean \
  --telescope MeerKAT \
  --observation-time 60 \
  --frequency-mhz 1300 \
  --pixels 512 \
  --catalog MIGHTEE \
  --imager wsclean \
  --clean-iterations 100 \
  --overwrite
```

Fast smoke run used for runtime verification:

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

Latest runtime result: completed with exit code 0. The manifest reports
`status: completed`, `catalog: MIGHTEE`, 9896 sky sources, completed simulation
and imaging milestones, visibility data, WSClean FITS/PNG image products, and
`weblog.html` under `smoke_mightee_wsclean_quick/`.

The earlier failed run used:

```text
conda run -n skasim skasim --telescope MeerKAT --seconds 60 --freq 1300 --pixels 1024 --catalog MIGHTEE --niter 500
```

It selected the default `oskar-dirty` imager, not WSClean, and failed during
the OSKAR simulation path. The output directory was
`20260523_064436_MeerKAT/`.

## Current CLI Access Layer

The 0.2 CLI now has one canonical spelling for each concept:

- sky model file: `--model`
- built-in catalog: `--catalog`
- generated Stokes I flux densities: `--flux-density`
- generated polarization: `--stokes-q`, `--stokes-u`, `--stokes-v`
- observing frequency: `--frequency-mhz`
- observation duration: `--observation-time`
- exact output directory: `--output-dir`
- WSClean iterations: `--clean-iterations`
- file-backed flux scaling: `--flux-scale`

Removed or renamed flags are hidden from `--help` and fail with targeted
migration messages: `--I`, `--stokes-i`, `--fits`, `--json`,
`--json-fg`, `--Q`, `--U`, `--V`, `--ref-freq`, `--freq`, `--bandwidth`,
`--delta-freq`, `--seconds`, `--prefix`, `--niter`, `--scale-I`,
`--cleaning`, and `--imaging-niter`.

JSON sky models are still supported through `--model path/to/sources.json`.
The separate `--json` and `--json-fg` entry points were removed because they
created parallel access paths for the same sky-source concept.

`--output-dir` is exact. Supplying `--output-dir smoke_mightee_wsclean` writes
to `smoke_mightee_wsclean/`; it does not append the telescope name.

## Review Notes

The development markdown is intentionally layered:

- Use this README as the collaborator summary.
- Use the PRD to understand the original rationale.
- Use the issue files to map implementation slices to acceptance criteria.
- Use the assessment and resolution files to understand why the strict
  post-assessment cleanup happened.
- Use `git log --oneline` and `git diff --name-only` for exact commit and file
  state in the working tree.
