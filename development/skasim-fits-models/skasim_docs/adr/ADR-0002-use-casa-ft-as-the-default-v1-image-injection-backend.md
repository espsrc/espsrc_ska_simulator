# ADR-0002: Use CASA `ft` as the default v1 image-injection backend

- Status: Accepted
- Date: 2026-05-30

## Context
The project needs a default backend for injecting image-based models into simulated visibilities. The backend must support:
- image-based prediction into a Measurement Set,
- additive multi-model injection,
- compatibility with the current requirement set for version 1,
- a path that does not depend on converting all image inputs into source lists.

WSClean predict was considered, but it is not suitable as the default backend because its standard prediction path is not aligned with the long-term requirements for smooth continuum treatment and additive multi-model injection.

DP3/WGridderPredict was also considered. It remains promising, especially for later performance-oriented work, but it was not chosen as the version-1 default.

## Decision
Use **CASA `ft`** as the default image-injection backend for version 1.

The workflow will:
- build a base Measurement Set,
- prepare backend-ready image inputs for CASA `ft`,
- inject one or more image-based models into `MODEL_DATA`,
- use additive injection semantics for multiple image-based entries,
- merge the result into final observed `DATA`.

## Rationale
CASA `ft` was chosen because it best matches the version-1 product requirements:
- It predicts from image-based models.
- It supports additive multi-model injection via `incremental=True`.
- It is a better fit than WSClean for the required composite-simulation workflow.
- It provides a stable, correctness-first default while the wider backend strategy remains open for later versions.

## Consequences
### Positive
- Supports additive multi-model injection in a straightforward way.
- Provides a single default backend for version 1.
- Aligns with the immediate need to combine catalog models and image-based models.

### Negative
- Introduces a CASA runtime dependency.
- May be slower than future alternatives on very large models.
- Requires helper conversion into CASA-ready model products.

## Alternatives considered
1. **WSClean predict as the default.** Rejected because it is not a sustainable default for the required image-injection semantics.
2. **DP3/WGridderPredict as the default.** Deferred. It remains a future backend candidate, particularly for performance-oriented development.
3. **Support multiple equivalent default backends from v1.** Rejected because it increases complexity without a clear product need in the first iteration.
