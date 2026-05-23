# Issue 04: Named Sky Model Sources

Type: AFK

## What to Build

Replace numeric catalog selection with named catalog selection and enforce one explicit sky model source per run. A sky model source can be file-backed, a named built-in catalog, or generated source intensities as fallback.

## Acceptance Criteria

- [ ] `--catalog MIGHTEE` selects the MIGHTEE built-in catalog.
- [ ] `--catalog GLEAM` selects the GLEAM built-in catalog.
- [ ] Numeric catalog values fail with a targeted migration message.
- [ ] Passing both a model file and a named catalog fails clearly.
- [ ] Generated source mode remains the fallback when no explicit sky model source is provided.
- [ ] Tests cover named catalogs, numeric migration failures, and mutually
  exclusive sky model source validation.

## Blocked By

- Issue 02: Pip-Only Install and Lazy Karabo Runtime

## Suggested Commit Message

`feat: use named sky model sources`
