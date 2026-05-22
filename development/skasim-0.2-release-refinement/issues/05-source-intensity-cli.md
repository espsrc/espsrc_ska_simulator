# Issue 05: Source Intensity CLI

Type: AFK

## What to Build

Replace unclear generated source intensity language around `--I` with clearer CLI flags. Generated source intensities should be expressed with `--flux-density` and `--stokes-i`, and they should only be valid in generated source mode.

## Acceptance Criteria

- [ ] `--flux-density` accepts one or more generated source intensities.
- [ ] `--stokes-i` is accepted as an alias for `--flux-density`.
- [ ] Source intensity flags create multiple generated sources when multiple values are provided.
- [ ] Source intensity flags are rejected when a file-backed or named catalogue sky model source is provided.
- [ ] `--I` is handled with the agreed 0.2 migration behavior.
- [ ] Tests cover generated source intensity behavior, aliases, and invalid combinations with explicit sky model sources.

## Blocked By

- Issue 04: Named Sky Model Sources

## Suggested Commit Message

`feat: clarify generated source intensity CLI`
