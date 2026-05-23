# Assessment 01 Resolution

This note records the implementation decisions made after reviewing
`assessment_01.md`.

## Accepted Findings

- Remove vestigial 0.1 configuration fields instead of preserving Python API
  compatibility.
- Make `source_flux_jy` the generated-source Stokes I source of truth.
- Reject `catalog=0` with the numeric catalog migration message.
- Remove stored `ImgConfig.algorithm`; `imager` is the image-product selector.
- Remove dead `args.I` and `cleaning` plumbing from CLI-to-config construction.
- Replace interactive visibility overwrite prompting with `FileExistsError`.
- Store failure milestone details as dictionaries.
- Render weblog duration cards from `simulation_*` and `imaging_*` milestones.
- Skip missing image files when rendering failed-run weblogs.
- Guard the sky-model summary template against missing `sky_model_loaded`.
- Use timezone-aware UTC datetimes in manifests.
- Use canonical WSClean `-auto-threshold`.
- Preserve run-scoped WSClean output files during temporary-file cleanup.
- Remove the stale commented FITS resolver block.
- Bump release metadata to `0.2.0`.
- Tighten tests for CLI/config migration, runtime imports, weblog failure
  rendering, WSClean cleanup, and release-doc path handling.

## Resulting Access Layer

- Public package exports are `ImgConfig`, `ObsConfig`, and `SimConfig`.
- `SimConfig` rejects extra fields, including removed 0.1 fields.
- `ObsConfig` and `ImgConfig` reject extra fields.
- Generated-source runs use `source_flux_jy`; default is `[10.0]`.
- Generated-source polarization remains available through `stokes_q_jy`,
  `stokes_u_jy`, and `stokes_v_jy`; supplied lists must match
  `source_flux_jy` length and default to zero when omitted.
- File-backed and named-catalog sky model runs have no generated source
  flux or polarization flags.
- CLI users use `--flux-density`, not `--I` or `--stokes-i`.
- CLI users use `--stokes-q`, `--stokes-u`, and `--stokes-v` for optional
  generated-source polarization, not uppercase `--Q`, `--U`, or `--V`.
- CLI users use `--imager wsclean`, not `--cleaning`.
- Built-in catalogs are named: `MIGHTEE`, `GLEAM`, and `SKAMid`.
- CLI users use `--model` for FITS, JSON, pickle, or Karabo model files; the
  separate `--fits`, `--json`, and `--json-fg` access paths were removed.
- CLI users use explicit unit-bearing options: `--frequency-mhz`,
  `--bandwidth-mhz`, `--channel-width-mhz`, and `--observation-time`.
- CLI users use `--output-dir` to name the exact output directory; the previous
  prefix-plus-telescope behavior was removed.
- CLI users use `--clean-iterations`, not `--niter`.

## Weblog And Preview Refinement

- APLpy and CMasher are accepted runtime dependencies for FITS image-product
  preview rendering.
- FITS previews are generated from in-memory 2D celestial HDUs so WCS axes,
  colorbar, beam overlay, and contours are handled by APLpy without temporary
  files.
- Preview display data are converted to mJy/beam and the rendered HDU metadata
  records `BUNIT=mJy/beam`.
- The weblog remains Jinja2-based. Removing Jinja2 in favor of large Python
  f-string HTML was rejected because it made the rendering code harder to
  maintain while adding heavier plotting dependencies anyway.
- The weblog now prioritizes science products: WSClean MFS Model, Clean, and
  Residual images are shown together, with the PSF exposed as a direct link.
  Dirty images are shown only when model/clean/residual products are absent.
- The telescope layout plot is displayed in the compact
  Observation/Telescope section instead of the supporting-plot gallery.
- Sky-model rendering now produces two explicit plot outputs, a full source
  model and a FoV-matched zoom, and the weblog displays both in the Sky Model
  section. These plots render source shapes as ellipses from major axis, minor
  axis, and position angle metadata, with a reversed colormap so bright sources
  are darker than faint sources. Compact sources below a FoV-dependent angular
  size threshold are rendered as flux-scaled crosses so the full-field plot
  remains complete, and source PA is converted from astronomical east-of-north
  convention to Matplotlib's x-axis convention.
- JSON sky-model loading now stores full source rows instead of reduced
  `(ra, dec, I)` rows so source sizes and position angles survive through
  `SkyModel.to_json()` and the sky-model preview renderer.
- The reference Gaussian catalog generator writes both JSON and DS9 FK5 region
  files. It now produces a broader demonstration population: most sources are
  centered around the 1 mJy scale, two deterministic bright sources are near
  60-100 mJy, several faint sources are in the 10-20 microJy range, and source
  sizes span sub-arcsec compact components through several-arcmin extended
  Gaussians. DS9 regions use ellipse definitions with the same source
  coordinates, major/minor axes, converted display PA, flux density, and
  spectral index metadata, plus cross markers so compact sources remain visible.
- Runtime timing is summarized in the small header metadata, and the detailed
  milestone timeline is placed after the science images.
- `--overwrite` now removes and recreates the exact output directory before
  logger/manifest setup. This prevents stale MeasurementSet lock files and
  WSClean products from previous runs from corrupting repeated smoke runs such
  as `demo_output/reference_meerkat`.
- CLI/pipeline execution forces Matplotlib's `Agg` backend and closes telescope
  plot figures after saving to prevent X11/ICE shutdown failures on headless
  systems.

## Verification

- `PYTHONPATH=/tmp/skasim-pydeps:src python -m py_compile src/skasim/*.py`
- `PYTHONPATH=/tmp/skasim-pydeps:src python -m pytest -q`
- Latest focused result:
  `109 passed in tests/test_config.py tests/test_cli.py tests/test_pipeline.py tests/test_manifest.py tests/test_imaging.py tests/test_runtime_imports.py tests/test_release_docs.py tests/test_weblog.py`.
- Runtime smoke result: the reference MeerKAT JSON-catalog WSClean command from
  `run_all.sh` completed with exit code 0 after the clean overwrite and
  non-interactive plotting fixes.
