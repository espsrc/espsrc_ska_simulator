# skasim 0.2 Implementation Assessment

Full review of the 12-commit implementation against the PRD and issue acceptance criteria.

---

## Executive Summary

The implementation delivers the planned features across all 12 issues. The architecture is sound — lazy Karabo imports, structured manifest, always-on weblog, named catalogs, explicit imager CLI. However, there are **significant dead code remnants**, **vestigial fields that should have been removed**, **test gaps**, and **several bugs**. Below is a per-issue breakdown followed by cross-cutting concerns.

---

## Issue 01: Development Baseline ✅ Clean

No code changes — docs-only commit. CONTEXT.md and development folder are well-structured.

---

## Issue 02: Pip-Only Install and Lazy Karabo Runtime

### [runtime.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/runtime.py)

**Good**: Clean, focused module. `require_karabo_module` is well-designed.

**Minor**: The conditional `if module_name.startswith("karabo") or exc.name == "karabo"` (line 25) silently re-raises non-Karabo ImportErrors. This is intentional but undocumented — add a comment.

### [sky.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/sky.py) — Fallback KaraboSkyModel

> [!WARNING]
> **Bug: Top-level imports break pip-only promise**

Lines 1-15 of `sky.py` unconditionally import `astropy.coordinates`, `astropy.table`, `radio_beam`, `xarray`, and `loguru` at module load time. While these are in `pyproject.toml` dependencies, `radio_beam` and `xarray` are heavy imports. More critically, **`sky.py` is imported by `pipeline.py`**, which is imported when running `skasim --help` if the CLI ever falls through to the `from .pipeline import run` path (line 221 of cli.py). The lazy import of pipeline only at execution time (line 221) is correct, but `__init__.py` imports from `manifest.py` which imports from `config.py` — this chain is safe. The `sky.py` chain is only pulled in at runtime. **Verdict: acceptable, but fragile.**

> [!IMPORTANT]  
> **Dead code: Fallback `KaraboSkyModel` at lines 470-486**

The fallback class duplicates a minimal Karabo API. It works, but `add_point_sources` uses `np.vstack` which will crash if `self.sources` has shape `(N, 3)` and new sources have shape `(M, 14)` — shape mismatch. This is never tested.

### Tests ([test_runtime_imports.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/tests/test_runtime_imports.py))

- `test_cli_help_does_not_require_karabo` — **Good**: subprocess isolation.
- `test_simulation_without_karabo_fails_with_installation_message` — **Good**: verifies the runtime error message.
- **Gap**: No test verifies that `import skasim` works without Karabo (only CLI help is tested).

---

## Issue 03: Conda Environment and Install Docs

### [environment.yml](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/environment.yml)

> [!WARNING]
> **Problem: `python=3.9` contradicts `pyproject.toml`'s `requires-python = ">=3.8.5"`**

Pick one truth. The conda env pins 3.9, which is fine for Karabo, but document why.

### [installation.rst](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/docs/installation.rst)

- Good separation of pip-only vs conda paths.
- WSClean Singularity documented as example, not default. ✅

### Tests ([test_release_docs.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/tests/test_release_docs.py))

- `test_environment_file_defines_skasim_conda_runtime` — verifies name, channels, karabo-pipeline. ✅
- `test_installation_docs_describe_pip_and_conda_paths` — asserts key phrases. ✅
- **Fragile**: These tests use relative paths (`Path("environment.yml")`) and will fail if tests are run from a different working directory. Should use `Path(__file__).resolve().parents[1]`.

---

## Issue 04: Named Sky Model Sources

### [config.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/config.py)

> [!CAUTION]
> **Bug: `_normalise_catalog` accepts `value == 0` as `None` (line 141)**

```python
if value in (None, "", 0):
    return None
```

This means `catalog=0` silently succeeds instead of triggering the numeric migration error. The PRD says *"Numeric catalog inputs should fail with targeted migration messages"*. A user passing `--catalog 0` would get silent fallback to generated sources instead of the migration message.

> [!WARNING]
> **Vestigial field: `cleaning: bool = False` still in SimConfig (line 136)**

This field was supposed to be *removed* in Issue 06. It's still present in the model, always hardcoded to `False` from the CLI (line 218). It should be deleted entirely — it's dead weight in every manifest JSON dump, confusing future readers.

> [!WARNING]
> **Redundant field: `algorithm` in ImgConfig (line 84)**

`algorithm` is always mechanically derived from `imager` (line 189 of cli.py: `"wsclean_clean" if args.imager == "wsclean" else "oskar_dirty"`). It carries no independent information. It should either be a computed property or removed, not a stored field that can get out of sync with `imager`.

### [cli.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/cli.py) — Catalog flag

- `--catalog` is the canonical built-in catalog flag. ✅
- **But**: `args.I` at line 203 refers to the deprecated `--I` flag. Since `_DeprecatedIAction` calls `parser.error()` (which exits), `args.I` is always `None` when execution reaches line 203. The fallback `args.I or [10.0]` is dead code that looks like it's still using the old flag. Should just be `[10.0]`.

### Tests

- Named catalogs, numeric migration, mutual exclusion — all well covered. ✅

---

## Issue 05: Source Intensity CLI

### [cli.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/cli.py)

- `--flux-density` and `--stokes-i` correctly map to `source_intensities`. ✅
- `--I` migration via `_DeprecatedIAction`. ✅

> [!IMPORTANT]
> **Problem: `source_intensities` vs `I` field confusion in pipeline**

In [pipeline.py line 245](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/pipeline.py#L245):
```python
source_intensities = config.source_intensities or config.I
```
This means the old `I` field is still the *actual* fallback for generated source mode. The `I` field default is `[10.0]` (config.py line 110). So if a user passes no intensity flags, `source_intensities` is `None`, and it falls back to `config.I = [10.0]`. This works but is confusing — the PRD wanted `--flux-density` to *replace* `--I`, yet the old field is still the source of truth for defaults. The `I` field should be removed from `SimConfig` and the default should live in `source_intensities` directly.

### Tests

- `test_flux_density_sets_generated_source_intensities` ✅
- `test_stokes_i_alias_sets_generated_source_intensities` ✅
- `test_legacy_i_flag_has_migration_message` ✅
- **Gap**: No test verifies that `source_intensities` with a catalog raises an error *at the CLI level* (only config-level test exists).

---

## Issue 06: Explicit Imager CLI

- `--imager oskar-dirty` default ✅
- `--imager wsclean` selection ✅
- `--cleaning` migration via `_DeprecatedCleaningAction` ✅

> [!NOTE]
> The `cleaning` field still exists in `SimConfig` (see Issue 04 above). This should have been removed in this issue.

### Tests — All good ✅

---

## Issue 07: Configurable WSClean Command

### [imaging.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/imaging.py)

- `build_wsclean_argv` uses `shlex.split` for the command. ✅
- `run_wsclean_command` uses `shell=False` and explicit `cwd`. ✅

> [!WARNING]
> **Bug: Mixed `--auto-threshold` vs `-auto-mask` flags (lines 103-106)**

```python
"--auto-threshold",   # double-dash
"0.3",
"-auto-mask",         # single-dash
"3",
```

WSClean uses single-dash for both: `-auto-threshold` and `-auto-mask`. The `--auto-threshold` with double-dash may work in some WSClean versions but is non-canonical and inconsistent.

> [!NOTE]
> **Hardcoded WSClean parameters** (lines 90-113)

`-multiscale`, `-mgain 0.8`, `-auto-mask 3`, `-channels-out 8`, `-join-channels`, `-local-rms` are all hardcoded. These should at minimum be documented as non-configurable defaults, and ideally some (like `niter`, which *is* configurable) should be consistently handled.

### Tests

- `test_build_wsclean_argv_uses_default_command` ✅
- `test_build_wsclean_argv_parses_singularity_command` ✅
- `test_run_wsclean_command_uses_argv_and_working_directory` ✅
- Good coverage of the core WSClean command contract.

---

## Issue 08: Stable WSClean Outputs

### [imaging.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/imaging.py)

- `wsclean_output_prefix` returns `{work_dir.name}_wsclean`. ✅
- `collect_wsclean_outputs` uses glob with prefix. ✅

> [!WARNING]
> **Bug: Old-file cleanup at lines 166-170 uses wrong glob pattern**

```python
for tmp in work_dir.glob("wsclean-00*.fits"):
```

This hardcoded pattern only cleans up files matching `wsclean-00*`. But WSClean temporary files could have different patterns. More importantly, this cleanup runs *before* `collect_wsclean_outputs`, so it could delete files that match the collection prefix if someone's prefix started with `wsclean-00`. This is fragile — cleanup should be explicit about what it removes.

### Tests

- `test_collect_wsclean_outputs_matches_only_configured_prefix` — **Good**: verifies old-run files are excluded. ✅
- `test_run_wsclean_imaging_uses_run_prefix_and_stable_outputs` — Integration-style test with mocks. ✅

---

## Issue 09: Structured Manifest Outputs

### [manifest.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/manifest.py)

- `OutputRecord` with kind, path, image_product_id, imager, role, metadata. ✅
- `RunManifest` with milestones, outputs, errors. ✅
- `model_dump_json` override for pretty printing. ✅

> [!WARNING]
> **`datetime.utcnow()` is deprecated in Python 3.12+**

Lines 71, 102 use `datetime.utcnow()`. This was deprecated in Python 3.12 in favor of `datetime.now(timezone.utc)`. While the project targets 3.8+, this will emit deprecation warnings on newer Pythons and should be updated.

> [!NOTE]
> **`model_dump_json` override hides Pydantic's native serialization**

The override at line 109 calls `json.dumps(self.model_dump(mode="json"), ...)` instead of using Pydantic's built-in `model_dump_json()`. This works but bypasses Pydantic's serialization optimizations and could cause subtle differences.

### Tests

- `test_manifest_serializes_structured_outputs` ✅
- `test_create_run_context_records_log_and_manifest_outputs` ✅
- `test_default_run_id_uses_second_precision` ✅

---

## Issue 10: Always-On Weblog

### [weblog.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/weblog.py)

> [!CAUTION]
> **Bug: `_find_image_outputs` crashes if image file doesn't exist (line 41)**

```python
data = fpath.read_bytes()  # FileNotFoundError if file missing
```

During a failed run, image files referenced in the manifest may not exist on disk. The function will crash with `FileNotFoundError`, preventing the failure weblog from rendering — defeating the purpose of "always-on weblog".

> [!WARNING]
> **Weblog phase duration uses wrong milestone names (lines 83-94)**

```python
if "phase_a_started" in milestone_lookup and "phase_a_completed" in milestone_lookup:
```

But the pipeline uses milestone names like `simulation_started`, `simulation_completed`, `imaging_started`, `imaging_completed` — not `phase_a_*`/`phase_b_*`. **These duration calculations will always be `None`** because the milestone names never match. This is dead code that silently produces no output.

### [weblog.html.j2](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/templates/weblog.html.j2)

> [!WARNING]
> **Jinja2 `first` filter at line 197 will error if no matching milestone exists**

```jinja2
{% set sky_ms = manifest.milestones | selectattr("name", "equalto", "sky_model_loaded") | first %}
```

Jinja2's `first` filter raises `StopIteration` (rendered as `UndefinedError`) if the sequence is empty. For failed runs where `sky_model_loaded` never happened, this crashes. Should use `| first | default(none)` or a conditional.

### Tests

- `test_weblog_renders_structured_outputs` — Only verifies output kind/path strings appear in HTML. **Very shallow**. Doesn't test:
  - Failed run weblog rendering
  - Image embedding
  - Milestone timeline rendering
  - Error section rendering

---

## Issue 11: Run Pipeline Records and Seconds Run IDs

### [pipeline.py](file:///mnt/scratch1/jmoldon/dev/espsrc_ska_simulator/src/skasim/pipeline.py)

> [!CAUTION]
> **Bug: `add_milestone` `details` parameter type mismatch (line 441)**

```python
ctx.add_milestone("simulation_failed", "failed", elapsed_s=time.time() - t_phase_a, details=str(exc))
```

`details` is typed as `Optional[dict]` in `Milestone`, but `str(exc)` is passed here. Pydantic will reject this — a string is not a dict. Same issue at line 459. This means **simulation/imaging failure milestones will crash** before the error can be recorded properly.

> [!WARNING]
> **`run_simulation` recomputes values already available from `build_observation`**

Lines 359-361 recompute `freq`, `fov`, and `delta_freq` from config, even though these were already computed in the `run()` orchestrator and passed to `build_observation`. This is the redundancy Issue 11 was supposed to fix.

> [!IMPORTANT]
> **`run_simulation` prompts for user input (line 353)**

```python
ans = input(f"{visibility_path} exists. Overwrite? (y/n): ")
```

Interactive `input()` in a library function is problematic. It blocks automated/scripted runs and tests. The `overwrite` flag should be the sole mechanism.

> [!NOTE]
> **Unused imports in pipeline.py**

- `glob` (line 5) — never used
- `pickle` (line 8) — used only in `_load_sky_from_file`
- `shutil` (line 9) — used only in `run_simulation`
- `sys` (line 10) — used only in `run_simulation`

These aren't wrong, but signal that pipeline.py is doing too much.

### Tests

- `test_run_uses_resolved_wsclean_imager` — **Good**: end-to-end pipeline test with mocks.
- `test_run_renders_weblog_on_failure` — **Good**: verifies failure path.
- `test_run_simulation_does_not_rebuild_observation` — **Good**: verifies no duplicate setup.
- **Gap**: No test for the `details=str(exc)` bug in failure milestones.

---

## Issue 12: Release Smoke Checks and Docs Polish

### Docs

- README uses 0.2 language ✅
- guide.rst comprehensive ✅
- examples.rst has named-catalog and FITS smoke shapes ✅

> [!WARNING]
> **pyproject.toml still says `version = "0.1.0"` (line 7)**

The PRD says 0.2 is the target. The version should be bumped to `0.2.0` or `0.2.0.dev0`.

### Tests

- `test_release_examples_use_0_2_cli_language` — **Clever**: asserts absence of deprecated flags. But `assert "--I" not in combined` will false-positive match against `--I` in words like `--Imager` (not currently an issue but fragile).

---

## Cross-Cutting Concerns

### 1. Dead / Vestigial Fields in SimConfig

| Field | Status | Action |
|---|---|---|
| `cleaning: bool = False` | Dead — always `False`, old `--cleaning` is deprecated | **Remove** |
| `I: List[float] = [10.0]` | Vestigial — superseded by `source_intensities` | **Remove**, move default to `source_intensities` |
| `algorithm` in ImgConfig | Redundant — derivable from `imager` | **Remove** or make `@property` |
| `sky_format` | Always `"auto"` from CLI (line 198) | Consider removing if never used |
| `source_names` | Never set anywhere, never read | **Remove** |

### 2. fits_helper.py — Large Commented-Out Block

Lines 24-228 contain a ~200-line commented-out block (`''' TODO: Review & simplify.`). This is dead code that should be either deleted or moved to a design doc. It includes `FitsColumnResolver`, `COLUMN_ALIASES`, etc. that are not used anywhere.

### 3. `__init__.py` Exports

```python
from .manifest import Milestone, OutputRecord, RunContext, RunManifest
```

This exports internal pipeline plumbing as the public API. `RunContext` especially is an implementation detail. Consider what the *public* API of skasim should be — likely just `SimConfig` and `run()`.

### 4. Test Quality Summary

| Test File | Coverage Quality | Issues |
|---|---|---|
| test_config.py | **Strong** | Good parametrized edge cases |
| test_cli.py | **Strong** | Covers all 0.2 CLI changes |
| test_imaging.py | **Strong** | Good mock-based argv testing |
| test_manifest.py | **Good** | Missing failure-path tests |
| test_pipeline.py | **Good** | Large, could be split; has the `details=str(exc)` gap |
| test_runtime_imports.py | **Adequate** | Missing bare `import skasim` test |
| test_sky.py | **Strong** | Thorough Source/SkyModel coverage |
| test_weblog.py | **Weak** | Single shallow test; no failure weblog test |
| test_release_docs.py | **Fragile** | Relies on CWD being repo root |
| test_utils.py | **Good** | Pre-existing, covers utilities well |

### 5. Consistency Issues

- **Milestone names**: Pipeline uses `simulation_started/completed/failed` and `imaging_started/completed/failed`, but weblog.py looks for `phase_a_started/completed` and `phase_b_started/completed`. These **never match**.
- **`details` type**: `Milestone.details` is `dict`, but pipeline passes `str(exc)` for failure milestones.
- **`show_exc` function**: Used in sky.py line 309 without import — will crash with `NameError` if an exception occurs during FITS table parsing.

---

## Priority Bug List

| # | Severity | Location | Description |
|---|---|---|---|
| 1 | 🔴 **High** | pipeline.py:441,459 | `details=str(exc)` — type mismatch crashes failure recording |
| 2 | 🔴 **High** | weblog.py:83-94 | Phase duration milestone names don't match pipeline names |
| 3 | 🟠 **Medium** | weblog.py:41 | `_find_image_outputs` crashes on missing files during failure weblog |
| 4 | 🟠 **Medium** | weblog.html.j2:197 | `| first` without default crashes on missing milestone |
| 5 | 🟠 **Medium** | config.py:141 | `catalog=0` silently accepted instead of migration error |
| 6 | 🟡 **Low** | sky.py:309 | `show_exc` not imported — NameError on FITS parse exception |
| 7 | 🟡 **Low** | imaging.py:103 | `--auto-threshold` double-dash inconsistent with WSClean |
| 8 | 🟡 **Low** | manifest.py:71,102 | `datetime.utcnow()` deprecated in Python 3.12+ |

## Priority Cleanup List

| # | Type | Location | Description |
|---|---|---|---|
| 1 | Dead code | config.py | Remove `cleaning`, `I`, `source_names`, `algorithm` fields |
| 2 | Dead code | fits_helper.py:24-228 | Remove 200-line commented-out block |
| 3 | Dead code | cli.py:203 | `args.I or [10.0]` — `args.I` is always None |
| 4 | Redundancy | pipeline.py:359-361 | `run_simulation` recomputes freq/fov/delta_freq |
| 5 | Design | pipeline.py:353 | Interactive `input()` in library code |
| 6 | Version | pyproject.toml:7 | Still says `0.1.0`, should be `0.2.0` |
| 7 | Tests | test_weblog.py | Only 1 shallow test — needs failure/milestone/error tests |
| 8 | Tests | test_release_docs.py | Uses relative paths, fragile CWD dependency |
