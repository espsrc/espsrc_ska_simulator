# ADR-0001: Separate model entry contracts from model preparation

- Status: Accepted
- Date: 2026-05-30

## Context
`skasim` already accepts component-style sky models through the Karabo/OSKAR path. The next stage of the project needs to support image-based inputs such as continuum maps, polarization maps, and later spectral cubes. These models often originate from different upstream tools and formats.

If the workflow tries to accept every upstream representation directly, the core pipeline becomes tightly coupled to many external formats and preparation conventions. That would make the product harder to reason about and harder to extend.

## Decision
The workflow will define a small number of **accepted model entry contracts**.

Model generation and most model preparation are outside the core pipeline. Helper scripts may be provided to convert common upstream products into the accepted contracts, but the workflow itself is responsible for:
- validating model entry type and required metadata,
- injecting the accepted model into the simulation,
- reporting provenance and assumptions.

For version 1, the accepted image-based contracts are:
- `continuum_i_alpha`
- `static_stokes_maps`

The existing component sky-model contract remains supported.

## Consequences
### Positive
- The workflow remains conceptually clean.
- New upstream model formats can be supported through helper conversion instead of core-pipeline redesign.
- The user-facing simulation interface stays stable.
- Future model families such as spectral cubes or dynamic models can be added as new contracts.

### Negative
- Some users will need a preparation step before running the workflow.
- Helper scripts become part of the practical user experience and need documentation.

## Alternatives considered
1. **Convert FITS images to source lists and treat them as ordinary Karabo sky models.** Rejected as the universal solution because it is not the right abstraction for polarization maps, spectral cubes, or future time-variable models.
2. **Accept any FITS-like model directly without typed contracts.** Rejected because it creates ambiguity and tightly couples the workflow to many heterogeneous formats.
