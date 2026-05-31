# Product Requirements Document

## Title
Image-based model injection for `skasim` v1

## Status
Draft for implementation planning

## Owner
Project owner: user / `skasim` maintainer

## Summary
`skasim` currently supports component-based sky models that are passed to Karabo/OSKAR and used to generate a Measurement Set. The next product increment is to support additional image-based model entry points without making model generation part of the core simulation workflow.

Version 1 adds support for:
- continuum image injection using `I(ν0)` + `alpha(x,y)`, and
- static polarization map injection using Stokes I/Q/U/V maps,

while preserving the existing component-sky-model path.

The workflow must also support composite simulations where more than one model contribution is combined in a single run.

---

## Problem statement
The current `skasim` pipeline assumes that the model passed to Karabo is catalog-like. This is insufficient for several realistic preparatory-science cases, including:
- injecting extended galaxies described by continuum images,
- injecting polarization maps,
- combining a background catalog with one or more image models,
- later extending to spectral cubes, time-variable models, and EoR-like signals.

A simple conversion from FITS images to source lists is not the desired long-term direction. It is not a good universal abstraction for continuum maps, polarization maps, spectral cubes, or future time-variable models. The pipeline therefore needs a distinct image-based injection path.

---

## Goals

### Primary goals
1. Support image-based model injection as a first-class workflow mode.
2. Keep model preparation outside the core simulation pipeline.
3. Preserve the current catalog/component sky-model workflow.
4. Allow multiple model contributions in one simulation run.
5. Produce a final observed Measurement Set as the main output.
6. Record relevant validation and provenance information in the weblog.

### Version 1 goals
1. Accept continuum `I+alpha` models.
2. Accept static Stokes I/Q/U/V map sets.
3. Support composite runs that combine:
   - an existing component sky model, and/or
   - one or more image-based models.
4. Use a single default image-injection backend that satisfies additive injection requirements.

---

## Non-goals for version 1
- Generating astrophysical models from first principles.
- Supporting arbitrary upstream image formats directly without normalization.
- Supporting spectral cubes in production.
- Supporting polarization cubes or RM/Faraday products.
- Supporting dynamic time-variable models.
- Supporting direction-dependent corruptions.
- Supporting mosaics.
- Providing multiple equivalent image-injection backends in v1.

---

## Success criteria
A version 1 implementation is successful if:
1. A user can provide a standard source-list sky model exactly as before.
2. A user can provide a continuum `I+alpha` image model and obtain a final MS containing the injected model plus noise.
3. A user can provide a static I/Q/U/V map set and obtain a final MS containing the injected model plus noise.
4. A user can combine a standard component sky model with one or more image models in the same run.
5. The workflow records what was injected, how it was interpreted, and what output files were produced.
6. Validation tests demonstrate correct additive behaviour, correct handling of model units, and correct handling of multi-model injection.

---

## Current state

### Existing execution flow
The current repository structure is:
- `src/skasim/cli.py`: CLI entry point.
- `src/skasim/config.py`: Pydantic validation.
- `src/skasim/pipeline.py`: orchestration.
- `src/skasim/sky.py` and `fits_helper.py`: sky-model loading.
- `src/skasim/imaging.py`: dirty/WSClean imaging.
- `manifest.py` and `weblog.py`: run manifest and report.

The current pipeline:
1. builds telescope and observation metadata,
2. builds a Karabo sky model,
3. runs Karabo/OSKAR simulation to produce a Measurement Set,
4. optionally performs imaging,
5. writes outputs and the weblog.

### Gap
The missing capability is a production-quality path for image-based model injection that is distinct from the catalog sky-model path.

---

## Users and use cases

### User types
- Astronomers simulating extended continuum sources.
- Astronomers simulating polarization maps.
- Users combining a wide-field source background with one or more special targets.
- Developers maintaining a general-purpose SKA/pathfinder simulation workflow.

### High-priority use cases
1. Inject a bright or extended galaxy from an `I+alpha` continuum model.
2. Combine a background source catalog with one image-based target.
3. Inject a static I/Q/U/V map set for polarization tests.
4. Run the full workflow and inspect the weblog for provenance and validation.

### Future use cases to keep compatible
1. Inject one or more spectral cubes.
2. Inject continuum plus line emission in the same run.
3. Inject time- and frequency-variable sources.
4. Inject EoR-like signals.

---

## Product decisions already made
1. Model generation is not part of the core pipeline.
2. The workflow accepts a limited set of standard model entry points.
3. The current catalog/component sky-model path remains supported.
4. Image-model units are standardized as Jy/pixel.
5. The default image-model assumption is intrinsic / PB-corrected sky.
6. The first version does not include direction-dependent effects.
7. The first version does not target mosaics.
8. The main required delivered data product is the final observed Measurement Set.
9. Static Stokes I/Q/U/V maps are in scope for version 1.
10. Spectral cubes, polarization cubes, RM products, dynamic models, and EoR-like signals are version 2+ items.
11. `catalog` is the canonical spelling in user-facing language and documentation.
12. The first implementation slice delivers `continuum_i_alpha` end-to-end and leaves `static_stokes_maps` schema-ready for the next phase.
13. Image-only runs first attempt a true empty/base MS and may fall back to an explicitly recorded zero-flux placeholder source if the runtime cannot create an empty MS.
14. `MODEL_DATA` is retained by default after final merge for validation and debugging.
15. FITS image-model previews must be shown in the weblog Sky Model area using the same preview rendering path as science image products.

---

## Core design decision

### Chosen image-injection backend for v1
Use **CASA `ft`** as the default image-based injection backend.

### Why this backend was chosen
- It supports image-based prediction into `MODEL_DATA`.
- It supports additive injection through `incremental=True`.
- It fits the need to combine multiple model contributions in one run.
- It is more compatible with the project requirements than WSClean as a default image-prediction path.

### Position of other backends
- **WSClean predict** is not the default backend because of spectral step/block behaviour and lack of suitable append semantics for the required long-term workflow.
- **DP3/WGridderPredict** remains a future backend candidate, primarily for performance-driven or later-phase development.

---

## Functional requirements

### FR1 — Existing component sky models must continue to work
The workflow must continue supporting the current catalog/component sky-model path with no regression in the existing user flow.

### FR2 — The workflow must support multiple model entries
The configuration must allow a run to contain zero or one component sky model and zero or more image-based model entries.

### FR3 — The workflow must support at least two image-based model entry types in v1
Supported image-based model entry types in v1:
1. `continuum_i_alpha`
2. `static_stokes_maps`

### FR4 — `continuum_i_alpha` must be an accepted user-facing model type
A `continuum_i_alpha` model must accept:
- a Stokes I map at a reference frequency,
- an `alpha` map on the same grid,
- the reference frequency value.

The workflow may internally convert this model into the format required by the injection backend.

### FR5 — `static_stokes_maps` must be an accepted user-facing model type
A `static_stokes_maps` model must accept Stokes I/Q/U/V maps defined on a common spatial grid. The initial implementation may allow a subset of these Stokes planes as long as the supplied set is explicit.

### FR6 — Sequential additive injection must be supported
The workflow must support additive multi-model injection. If more than one image-like model is provided, the predicted visibilities must be added rather than replaced.

### FR7 — Catalog models and image models must be combinable
A simulation run must be able to combine:
- the existing Karabo component sky model path, and
- one or more image-based model entries.

### FR8 — Final observed `DATA` must be produced
The delivered Measurement Set must contain final observed `DATA`. Retaining `MODEL_DATA` is optional but recommended for debugging and validation.

### FR9 — Validation and provenance must be recorded
The weblog and/or manifest must record:
- input model type(s),
- input file paths,
- unit interpretation,
- reference frequency where relevant,
- Stokes content,
- additive injection order,
- final output products.

### FR10 — Clear failure modes must exist
The workflow must fail early and clearly when:
- required model files are missing,
- image dimensions or WCS do not match expected assumptions for a model entry,
- model units are incompatible with the accepted contract,
- required metadata such as reference frequency is missing for `continuum_i_alpha`.

---

## Model entry contracts

### Contract A — `component_sky_model`
Purpose: existing path.

Accepted content:
- source list / catalog / Gaussian-component style input already supported by the current pipeline.

Injection path:
- Karabo/OSKAR directly.

### Contract B — `continuum_i_alpha`
Purpose: continuum extended-source injection.

Required content:
- Stokes I FITS image in Jy/pixel,
- `alpha` FITS image on the same grid,
- reference frequency.

Expected semantics:
- the continuum spectrum is defined as a per-pixel power law around the reference frequency.

Version 1 backend path:
- helper conversion into the form required by CASA `ft`, then prediction into `MODEL_DATA`.

### Contract C — `static_stokes_maps`
Purpose: static polarization-map injection.

Required content:
- one or more of Stokes I, Q, U, V maps,
- common spatial grid,
- Jy/pixel units.

Expected semantics:
- static spatial polarization structure for version 1.

Version 1 backend path:
- helper conversion into the form required by CASA `ft`, then prediction into `MODEL_DATA`.

---

## Proposed workflow changes

### High-level orchestration
1. Parse configuration.
2. Validate model entry types and files.
3. Build telescope and observation metadata.
4. If a component sky model is supplied, generate the initial MS using Karabo/OSKAR.
   - If no component sky model is supplied, generate an otherwise valid base MS for later image injection.
5. Initialize or clear the model column used for image-based injection.
6. For each image-based model entry in the configured order:
   - prepare the backend-ready form,
   - inject into `MODEL_DATA` using CASA `ft`,
   - use additive mode for all but possibly the first image-based injection.
7. Merge model visibilities into `DATA` to produce final observed data.
8. Run requested imaging.
9. Write manifest and weblog.

### Order of operations
The order of multiple image-based injections must be deterministic and recorded in the run report.

---

## Configuration requirements

### New top-level concept
The configuration must support **multiple model entries**.

A representative shape is:

```json
{
  "models": [
    {
      "type": "component_sky_model",
      "path": "background_catalog.fits"
    },
    {
      "type": "continuum_i_alpha",
      "stokes_i": "galaxy_i0.fits",
      "alpha": "galaxy_alpha.fits",
      "reference_frequency_hz": 1.4e9
    },
    {
      "type": "static_stokes_maps",
      "stokes_i": "pol_I.fits",
      "stokes_q": "pol_Q.fits",
      "stokes_u": "pol_U.fits",
      "stokes_v": "pol_V.fits"
    }
  ]
}
```

This example is illustrative. Exact field names may differ, but the design must support multiple typed model entries.

### Config validation requirements
The Pydantic config layer must validate:
- model entry type,
- required fields per type,
- file existence,
- mutually dependent fields,
- forbidden fields for a type,
- ordering preservation.

---

## Pipeline/module impact

### `config.py`
Add schema support for a list of typed model entries.

### `pipeline.py`
Add orchestration for:
- building a base MS,
- preparing image-model inputs,
- invoking CASA `ft`,
- additive sequential injection,
- final merge into `DATA`.

### `sky.py` / `fits_helper.py`
Keep the existing component-sky-model path intact. Add validation utilities or shared helpers only if they do not blur the distinction between catalog loading and image-model preparation.

### New model-preparation helpers
Add helper functionality to normalize accepted image-model entry types into backend-ready inputs. These may live in a dedicated module rather than inside the existing catalog loader.

### `imaging.py`
No fundamental redesign is required. Imaging remains downstream of visibility generation.

### `manifest.py` / `weblog.py`
Extend run reporting to capture model-entry metadata and injection provenance.

---

## Reporting requirements
The weblog should include a model-injection section with:
- the list of model entries,
- the declared model type of each entry,
- the injection order,
- any intermediate files created for injection,
- the backend used,
- the final MS path,
- warnings or assumptions applied.
- FITS image-model previews in the current Sky Model area, rendered with the same FITS preview configuration used for science image previews.

Recommended additional reporting:
- model file dimensions,
- Stokes present,
- reference frequency,
- unit checks,
- additive merge summary.

---

## Validation requirements

### Required validation tests
1. **Regression test:** existing component sky-model workflow still works.
2. **Single `I+alpha` test:** inject one continuum image model and confirm a valid final MS is produced.
3. **Single Stokes-map test:** inject one static I/Q/U/V model set and confirm a valid final MS is produced.
4. **Composite test:** component sky model + one `I+alpha` image model.
5. **Multi-image additive test:** two image-based model entries injected sequentially.
6. **Ordering test:** prove that the configured injection order is preserved and reported.
7. **Metadata test:** weblog and manifest record the expected model-entry provenance.

### Recommended correctness checks
1. Confirm that additive multi-model injection changes the final visibilities relative to single-model injection.
2. Confirm that missing required fields cause clear failures.
3. Confirm that unsupported model types fail with a clear message.
4. Confirm that unit assumptions are enforced consistently.

---

## Risks and constraints

### R1 — CASA dependency
The solution depends on CASA `ft`. Packaging, invocation, and environment stability must be handled carefully.

### R2 — Performance for large models
CASA `ft` may be slower than later alternatives for large images. This is acceptable for v1 if functionality and correctness are achieved.

### R3 — Distinction between workflow and helper scripts
It must remain clear which tasks are core workflow responsibilities and which are preparation/helper responsibilities.

### R4 — Version-1 polarization scope
Static I/Q/U/V maps are in scope, but spectral-polarization behaviour is not. Documentation must be explicit to avoid confusion.

### R5 — Base-MS creation without a catalog model
If the user supplies only image-based models, the workflow still needs a valid base MS. The implementation must define how that empty/base MS is generated.

---

## Resolved implementation decisions
1. The canonical configuration shape is a typed `models` list.
2. Image-model preparation helpers live outside the catalog/component sky-model loader.
3. Backend-ready CASA files may be persisted in the run directory when useful for provenance and debugging.
4. `MODEL_DATA` is kept by default after final merge.
5. Static Stokes map backend representation is deferred until the phase-2 implementation.

---

## Out-of-scope but important follow-ons
1. Spectral cubes.
2. Continuum + line composite modes with cube injection.
3. Polarization cubes.
4. RM/Faraday products.
5. Time-variable models.
6. EoR signal injection.
7. DP3 backend evaluation and benchmarking.

---

## Rollout recommendation

### Phase 1
- introduce typed model entries,
- support `continuum_i_alpha`,
- support additive CASA `ft` injection,
- support composite catalog + image runs,
- update weblog/manifest.

### Phase 2
- add `static_stokes_maps`,
- expand validation/reporting,
- refine tests and examples.

### Phase 3
- plan v2 support for spectral cubes and more advanced model classes.

---

## Acceptance checklist
- [ ] Existing catalog workflow preserved.
- [ ] Multiple model entries supported.
- [ ] `continuum_i_alpha` supported.
- [ ] `static_stokes_maps` supported.
- [ ] Sequential additive injection supported.
- [ ] Final observed `DATA` produced.
- [ ] Weblog/manifest updated.
- [ ] Regression and validation tests added.
