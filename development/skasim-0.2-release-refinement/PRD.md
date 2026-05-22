# PRD: skasim 0.2 Release Refinement and Architecture Cleanup

## Problem Statement

skasim 0.1 works as an early release baseline, but the repository is not yet shaped for a robust 0.2 release. Users need a clearer command-line interface, reliable installation paths, configurable WSClean execution, inspectable run records, and cleaner module structure without losing the current scientific functionality.

The current interface still exposes ambiguous or transitional concepts: numeric catalogue IDs, `--I` as the source intensity flag, `--cleaning` as a boolean imager switch, and a hardcoded `wsclean` invocation. The code also imports Karabo in places that prevent pip-only imports and lightweight checks from working. Run records are partially flat, weblog generation only happens on success, and imaging output naming/discovery is brittle for future multi-image-product runs.

The 0.2 release should keep the useful 0.1 behavior, but make the simulator easier to install, test, run, inspect, and extend.

## Solution

Create a 0.2 release refinement that preserves the core run behavior while replacing ambiguous interfaces with explicit domain language.

A run should accept one sky model source: a file-backed sky model, a named built-in catalogue, or generated source intensities. Built-in catalogues should be selected by name rather than numeric ID. Generated source mode should use clearer source intensity flags. Imaging should be selected by an explicit imager name, with WSClean configured through a WSClean command that defaults to `wsclean` but can be set to a Singularity invocation on this machine.

The package should support pip-only installation for imports, configuration validation, CLI help, docs, and lightweight tests. Full simulation runs remain supported through a conda environment that installs Karabo. Runtime dependencies on Karabo should fail only when execution actually needs Karabo, with clear messages.

Every run should produce a manifest and weblog, including failed runs. The manifest should record structured output records so image products can be grouped by kind, imager, role, and path. The implementation should deepen the sky model resolution, imaging, runtime backend, run pipeline, and configuration language modules enough to make behavior testable through stable interfaces.

## User Stories

1. As a simulator user, I want to install skasim with pip for lightweight use, so that I can inspect the CLI and validate configuration without a full Karabo environment.
2. As a simulator user, I want a conda environment file for full simulations, so that I can reproduce the supported Karabo installation path.
3. As a simulator user, I want installation instructions for both pip-only and conda workflows, so that I can choose the right setup for my task.
4. As a simulator user, I want `skasim --help` to work without Karabo installed, so that I can learn the interface before setting up the full runtime stack.
5. As a simulator user, I want clear runtime errors when Karabo is missing, so that I know how to fix my environment.
6. As a simulator user, I want built-in catalogues selected by name, so that commands and logs are readable.
7. As a simulator user, I want numeric catalogue IDs removed, so that outdated commands do not silently encode unclear behavior.
8. As a simulator user, I want old numeric catalogue IDs to fail with a migration message, so that I know which named catalogue to use.
9. As a simulator user, I want `--catalog` accepted as an alias for `--catalogue`, so that either spelling works at the CLI.
10. As a simulator user, I want `--catalogue MIGHTEE`, so that I can use Karabo's built-in MIGHTEE catalogue source.
11. As a simulator user, I want `--catalogue GLEAM`, so that I can use Karabo's built-in GLEAM catalogue source.
12. As a simulator user, I want file-backed sky models to remain available through `--model`, so that I can run from FITS, JSON, pickle, or Karabo model files.
13. As a simulator user, I want `--model MIGHTEE_DR1.fits` to be treated differently from `--catalogue MIGHTEE`, so that file-backed and built-in catalogue sources remain distinct.
14. As a simulator user, I want exactly one explicit sky model source in 0.2, so that accidentally passing both `--model` and `--catalogue` fails rather than being silently prioritized.
15. As a future simulator user, I want the internal design to be ready for multiple sky model sources, so that later releases can combine background catalogues and foreground source models.
16. As a simulator user, I want generated source intensity mode to remain available when no file or catalogue is supplied, so that quick random-source runs still work.
17. As a simulator user, I want a clearer source intensity flag than `--I`, so that I understand it creates generated sources rather than an I/Q/U/V vector.
18. As a simulator user, I want `--flux-density` to define generated source intensities, so that random source commands read naturally.
19. As a simulator user, I want `--stokes-i` as a clear alias for source intensities, so that the Stokes meaning is explicit.
20. As a simulator user, I want old `--I` usage handled with a migration path, so that existing commands fail clearly or remain intentionally supported only if decided during implementation.
21. As a simulator user, I want source intensity flags rejected when a file or catalogue sky model source is provided, so that unused CLI values do not hide mistakes.
22. As a simulator user, I want an explicit `--imager` option, so that image production is selected by imager name rather than a cleaning boolean.
23. As a simulator user, I want `--imager oskar-dirty` as the default, so that existing dirty-image behavior remains the default run mode.
24. As a simulator user, I want `--imager wsclean`, so that WSClean cleaned imaging is selected explicitly.
25. As a simulator user, I want `--cleaning` removed with a migration message, so that old commands point me to `--imager wsclean`.
26. As a simulator user, I want to configure the WSClean command, so that this machine can run WSClean through Singularity.
27. As a simulator user on this machine, I want to pass `--wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"`, so that I can use the available WSClean container.
28. As a simulator user, I want `wsclean` to remain the default WSClean command, so that normal installations do not require extra configuration.
29. As a simulator maintainer, I want WSClean command execution to avoid shell-specific behavior, so that command construction is safer and easier to test.
30. As a simulator maintainer, I want WSClean output discovery to use a stable run-specific prefix, so that old files in the work directory are not mistaken for current outputs.
31. As a simulator user, I want output filenames to be stable and run-scoped, so that I can predict where outputs are written.
32. As a simulator user, I want scientific and configuration metadata recorded in the manifest rather than encoded in filenames, so that output names do not become brittle.
33. As a simulator user, I want every run to produce a weblog, so that results are always inspectable.
34. As a simulator user, I want failed runs to produce a weblog, so that I can inspect milestones and errors without reading raw logs first.
35. As a simulator user, I do not want a disable flag for the weblog, so that every run has a consistent human-readable record.
36. As a simulator user, I want the manifest to record structured outputs, so that visibility data, image products, plots, logs, and weblogs can be distinguished.
37. As a future simulator user, I want the manifest to support multiple image products, so that later releases can create several images with different imaging configurations from one run.
38. As a simulator user, I want 0.2 to produce one image product through the CLI, so that the release remains focused.
39. As a simulator maintainer, I want the internal imaging design to accept a list of image product requests, so that multiple image products can be added later without rewriting image production.
40. As a simulator user, I want default run IDs to include seconds, so that repeated runs in the same minute do not collide.
41. As a simulator user, I do not want random suffixes in default run IDs, so that output directories remain readable.
42. As a simulator user, I want existing `--output-prefix` behavior to remain understandable, so that I can still choose a run directory base name.
43. As a simulator user, I want CLI flags for random source mode, FITS catalogue ingestion, WSClean imaging, and weblog output to be documented, so that release examples are executable.
44. As a simulator maintainer, I want the sky model resolver tested through its interface, so that file, catalogue, and generated-source behavior can change internally without breaking callers.
45. As a simulator maintainer, I want Karabo runtime imports isolated, so that lightweight tests can run without the full simulation stack.
46. As a simulator maintainer, I want image production tested without executing WSClean, so that command construction and output collection are reliable.
47. As a simulator maintainer, I want run pipeline tests to verify success and failure records, so that manifest and weblog behavior is stable.
48. As a simulator maintainer, I want configuration and CLI tests to cover migration messages, so that breaking 0.2 changes are intentional and helpful.
49. As a simulator maintainer, I want installation docs to match the checked-in environment file, so that release setup does not drift.
50. As a simulator maintainer, I want the architecture to avoid unnecessary broad rewrites, so that 0.2 improves depth where it supports concrete release needs.

## Implementation Decisions

- Treat version 0.1 as the current baseline and this PRD as the 0.2 release refinement.
- Keep the current scientific functionality: random source mode, FITS catalogue ingestion, WSClean imaging, and weblog output.
- Preserve file-backed sky model behavior while clarifying that a file-backed sky model and a built-in catalogue are different sky model source types.
- Build or modify a deep sky model resolver module. Its interface should resolve one sky model and phase center from one sky model source. Internally it should be shaped so later releases can compose multiple sky model sources.
- Keep `skasim`'s current `SkyModel` inheritance approach for 0.2 unless implementation reveals a concrete blocker.
- Remove numeric built-in catalogue IDs from the 0.2 interface.
- Built-in catalogues must be selected by name, currently `MIGHTEE`, `GLEAM`, and, if retained, `SKAMid`.
- Numeric catalogue inputs should fail with targeted migration messages rather than silently mapping to names.
- Keep `catalogue` as the canonical domain term. Allow `catalog` only as a CLI spelling alias.
- Enforce one explicit sky model source per run in 0.2. Passing both a model file and a named catalogue should be invalid.
- Treat generated source intensity mode as the fallback when no file or catalogue sky model source is provided.
- Reject generated source intensity flags when a file or catalogue sky model source is provided.
- Replace unclear source intensity language around `--I` with clearer `--flux-density` and `--stokes-i` flags.
- Do not overload a single flux flag to mean either multiple generated source intensities or an I/Q/U/V vector.
- Build or modify a runtime backend module that isolates Karabo imports behind execution-time seams.
- Pip-only installation should support imports, configuration validation, CLI help, docs, and lightweight tests.
- Full simulations remain supported through conda with Karabo installed.
- Add an `environment.yml` for the supported full simulation environment named `skasim`.
- Do not add a pip extra for Karabo unless it is verified to work; prefer the conda environment as the supported full-runtime path.
- Update installation documentation to include both pip-only and conda workflows.
- Build or modify an image production module that treats image products as a list internally, while the 0.2 CLI still creates one image product.
- Replace `--cleaning` with an explicit `--imager` option.
- Default imager should be `oskar-dirty`.
- WSClean imaging should be selected through `--imager wsclean`.
- Removed `--cleaning` usage should fail with a targeted migration message pointing to `--imager wsclean`.
- Add a WSClean command configuration field and CLI flag.
- The WSClean command should default to `wsclean`.
- The WSClean command must support this machine's Singularity invocation.
- WSClean execution should parse the configured command into argv and run without shell execution.
- WSClean execution should use an explicit working directory rather than changing process-global current directory.
- WSClean output should use a stable run-specific `-name` prefix.
- WSClean output collection should collect only files matching the run-specific prefix.
- Move scientific/configuration meaning out of output filenames and into the manifest.
- Keep output filenames stable and predictable.
- Build or modify the run pipeline module so a run always saves manifest state and renders the weblog, including on failure.
- Weblog is always on and cannot be disabled.
- Add structured output records to the manifest.
- Structured output records should distinguish at least visibility data, image products, plots, logs, manifest, and weblog.
- Structured image product records should include enough metadata to group by image product and imager.
- Default run IDs should use second precision: `YYYYMMDD_HHMMSS_<telescope>`.
- Do not add random suffixes to run IDs.
- Keep broad refactors limited to changes directly justified by 0.2 release needs.

## Testing Decisions

- Tests should focus on external behavior at stable module interfaces rather than implementation details.
- Configuration and CLI tests should verify named catalogue behavior, `--catalog` aliasing, numeric catalogue migration failures, imager selection, removed `--cleaning` migration, source intensity flags, and WSClean command configuration.
- Sky model resolver tests should cover file-backed, named catalogue, and generated source intensity modes. Karabo catalogue behavior should be isolated behind test adapters where possible.
- Runtime backend tests should verify that pip-only imports and CLI help do not require Karabo, and that execution without Karabo fails with a clear message.
- Image production tests should verify WSClean argv construction, working-directory handling, stable output prefix use, and output collection without running WSClean.
- Run pipeline tests should verify manifest persistence, structured output records, success weblog generation, and failed-run weblog generation.
- Manifest/weblog tests should verify that structured outputs render correctly and that failed runs show error information.
- Installation checks should verify that the pip-only path supports import/help/lightweight tests, and that the conda environment file is documented as the full simulation path.
- Existing tests in the repository provide prior art for config validation, source and sky model behavior, pipeline helpers, and utility behavior.
- Smoke testing should include a MeerKAT run with a named MIGHTEE catalogue and WSClean imager using the configured WSClean command where the environment supports it.

## Out of Scope

- Dockerfile support is out of scope for this PRD and will be handled later.
- Multiple sky model source composition is out of scope for 0.2, although the internal resolver should not block it later.
- Multiple image products exposed through the CLI are out of scope for 0.2, although the internal imaging design should be ready for it later.
- Replacing skasim's current `SkyModel` representation with a non-Karabo representation is out of scope unless implementation reveals a concrete blocker.
- Broad scientific default changes are out of scope unless required by a bug fix.
- Replying to or resolving GitHub PR review comments is out of scope because there is no existing PR for this work.
- Docker or containerized full-stack release automation is out of scope.

## Further Notes

- The repository currently has no `environment.yml`; one should be added for 0.2.
- The current base environment could not run tests because dependencies such as Pydantic were missing and the package was not installed.
- Syntax compilation of current Python files passed before implementation began.
- `gh` is not installed in the local environment, but GitHub connector tooling is available for issue creation.
- The current WSClean container available on this machine is `/mnt/software/containers/wsclean-3.10-dysco.sif`.
- The architecture review generated during planning is available locally at `/tmp/architecture-review-20260522-154154.html`.
- Domain language captured in `CONTEXT.md` should be used in implementation and documentation: catalogue, sky model, sky model source, source intensity, run, image product, WSClean command, and weblog.
