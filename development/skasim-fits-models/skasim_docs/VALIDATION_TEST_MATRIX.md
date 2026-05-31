# Validation Test Matrix

## Objective
Define the minimum validation set for image-based model injection in `skasim` version 1.

---

## A. Regression and workflow integrity

### A1. Existing component-sky-model regression
- **Input**: standard catalog/component model only.
- **Expected result**: existing workflow still produces a valid MS and downstream outputs.
- **Purpose**: protect current behaviour.

### A2. Image-only workflow bootstrapping
- **Input**: image-based model only, no component catalog.
- **Expected result**: workflow produces a valid base MS, injects the image model, and writes final observed `DATA`.
- **Purpose**: validate base-MS creation without a source catalog.

---

## B. Continuum `I+alpha`

### B1. Single `I+alpha` injection
- **Input**: one `continuum_i_alpha` model.
- **Expected result**: successful run and valid final MS.
- **Purpose**: core functionality.

### B2. Missing reference frequency failure
- **Input**: one `continuum_i_alpha` model without reference frequency.
- **Expected result**: clear validation failure.
- **Purpose**: metadata enforcement.

### B3. WCS mismatch failure
- **Input**: Stokes I and alpha maps on mismatched grids.
- **Expected result**: clear validation failure.
- **Purpose**: protect against invalid model pairing.

---

## C. Static Stokes maps

### C1. Single static Stokes-map set
- **Input**: static Stokes I/Q/U/V maps.
- **Expected result**: successful run and valid final MS.
- **Purpose**: version-1 polarization scope.

### C2. Partial Stokes set
- **Input**: only I/Q/U or only I.
- **Expected result**: accepted if the chosen contract permits subsets and the weblog records what was supplied.
- **Purpose**: verify subset handling.

### C3. Stokes-grid mismatch failure
- **Input**: Stokes maps with differing dimensions or WCS.
- **Expected result**: clear validation failure.
- **Purpose**: enforce contract integrity.

---

## C2. CASA Taylor terms

### C2.1. Single CASA Taylor-term pair
- **Input**: one `casa_taylor_terms` model with `tt0` and `tt1`.
- **Expected result**: workflow passes the existing CASA image tables directly to CASA `ft`, merges `MODEL_DATA` into `DATA`, and records the model paths.
- **Purpose**: cover prepared CASA image fixtures such as local `*.model.tt0/.tt1` products.

### C2.2. Missing CASA table content
- **Input**: `casa_taylor_terms` path that exists but is not a CASA image table.
- **Expected result**: clear validation failure.
- **Purpose**: avoid passing arbitrary directories into CASA `ft`.

---

## D. Composite and additive injection

### D1. Component catalog + `I+alpha`
- **Input**: one component sky model and one `continuum_i_alpha` model.
- **Expected result**: successful composite simulation.
- **Purpose**: mixed-mode support.

### D2. Two image-based model entries
- **Input**: two image-based models injected sequentially.
- **Expected result**: both contributions are present in the final data.
- **Purpose**: additive injection.

### D3. Ordered multi-model injection
- **Input**: multiple image-based models in a known configured order.
- **Expected result**: injection order is preserved and reported.
- **Purpose**: deterministic orchestration.

---

## E. Provenance and reporting

### E1. Weblog content
- **Input**: representative composite run.
- **Expected result**: weblog lists model entries, types, order, backend, and output paths.
- **Purpose**: transparency and traceability.

### E2. Manifest content
- **Input**: representative composite run.
- **Expected result**: manifest captures generated files and key model-injection metadata.
- **Purpose**: machine-readable provenance.

### E3. FITS model preview
- **Input**: one `continuum_i_alpha` model.
- **Expected result**: weblog Sky Model section displays a FITS-model preview rendered through the same FITS preview path used for science image products.
- **Purpose**: make injected image models visually inspectable in the run report.

---

## F. Negative tests

### F1. Unknown model type
- **Input**: unsupported `model.type`.
- **Expected result**: clear error.

### F2. Missing model file
- **Input**: path to a non-existent FITS file.
- **Expected result**: clear error.

### F3. Invalid unit handling
- **Input**: model image with incompatible or missing unit metadata.
- **Expected result**: behaviour consistent with the implementation decision; either clear failure or explicit, recorded assumption.

---

## Recommended synthetic fixtures
- small catalog-only test field,
- small continuum Gaussian image + alpha map,
- small static polarization map set,
- composite test combining the catalog and a continuum image model.
