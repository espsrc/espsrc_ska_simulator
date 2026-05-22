# Issue 04: Named Sky Model Sources

Type: AFK

## What to Build

Replace numeric catalogue selection with named catalogue selection and enforce one explicit sky model source per run. A sky model source can be file-backed, a named built-in catalogue, or generated source intensities as fallback.

## Acceptance Criteria

- [ ] `--catalogue MIGHTEE` selects the MIGHTEE built-in catalogue.
- [ ] `--catalogue GLEAM` selects the GLEAM built-in catalogue.
- [ ] `--catalog` is accepted as a CLI alias for `--catalogue`.
- [ ] Numeric catalogue values fail with a targeted migration message.
- [ ] Passing both a model file and a named catalogue fails clearly.
- [ ] Generated source mode remains the fallback when no explicit sky model source is provided.
- [ ] Tests cover named catalogues, numeric migration failures, alias behavior, and mutually exclusive sky model source validation.

## Blocked By

- Issue 02: Pip-Only Install and Lazy Karabo Runtime

## Suggested Commit Message

`feat: use named sky model sources`
