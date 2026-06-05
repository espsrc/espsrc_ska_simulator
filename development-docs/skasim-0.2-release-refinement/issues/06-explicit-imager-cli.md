# Issue 06: Explicit Imager CLI

Type: AFK

## What to Build

Replace the `--cleaning` boolean with an explicit imager choice. The CLI should express image product creation through `--imager oskar-dirty` or `--imager wsclean`.

## Acceptance Criteria

- [ ] `--imager oskar-dirty` is the default.
- [ ] `--imager wsclean` selects WSClean imaging.
- [ ] `--cleaning` fails with a targeted migration message pointing to `--imager wsclean`.
- [ ] The resolved imager is recorded in configuration and run records.
- [ ] Tests cover default imager behavior, WSClean imager behavior, and removed `--cleaning` migration.

## Blocked By

- Issue 02: Pip-Only Install and Lazy Karabo Runtime

## Suggested Commit Message

`feat: replace cleaning flag with explicit imager`
