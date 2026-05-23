# Issue 05: Source Intensity CLI

Type: AFK

## What to Build

Replace unclear generated source intensity language around `--I` with clearer
CLI flags. Generated source Stokes I flux densities should be expressed with
`--flux-density`; optional generated-source polarization should be expressed
with `--stokes-q`, `--stokes-u`, and `--stokes-v`. These flags are only valid
in generated source mode.

## Acceptance Criteria

- [ ] `--flux-density` accepts one or more generated source Stokes I flux
  densities.
- [ ] `--stokes-i` fails with a targeted migration message for
  `--flux-density`; it is not a public alias after the post-assessment CLI
  consolidation.
- [ ] `--stokes-q`, `--stokes-u`, and `--stokes-v` accept optional generated
  source polarization values and require the same length as `--flux-density`.
- [ ] Source flux flags create multiple generated sources when multiple values are provided.
- [ ] Source flux and polarization flags are rejected when a file-backed or named catalog sky model source is provided.
- [ ] `--I` is handled with the agreed 0.2 migration behavior.
- [ ] Tests cover generated source flux behavior, removed aliases,
  polarization length validation, and invalid combinations with explicit sky
  model sources.

## Blocked By

- Issue 04: Named Sky Model Sources

## Suggested Commit Message

`feat: clarify generated source intensity CLI`
