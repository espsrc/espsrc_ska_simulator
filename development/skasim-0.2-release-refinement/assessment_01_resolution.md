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
- The telescope layout plot is displayed in the compact
  Observation/Telescope section instead of the supporting-plot gallery.
- Runtime timing is summarized in the small header metadata, and the detailed
  milestone timeline is placed after the science images.
- The broader pipeline-level sky-model preview generator proposed during
  weblog polishing was rejected for now; sky-model preview generation needs a
  separate design and tests before it becomes pipeline behavior.

## Verification

- `PYTHONPATH=/tmp/skasim-pydeps:src python -m py_compile src/skasim/*.py`
- `PYTHONPATH=/tmp/skasim-pydeps:src python -m pytest -q`
- Latest focused result:
  `109 passed in tests/test_config.py tests/test_cli.py tests/test_pipeline.py tests/test_manifest.py tests/test_imaging.py tests/test_runtime_imports.py tests/test_release_docs.py tests/test_weblog.py`.
