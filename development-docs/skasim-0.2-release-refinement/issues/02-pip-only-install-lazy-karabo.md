# Issue 02: Pip-Only Install and Lazy Karabo Runtime

Type: AFK

## What to Build

Allow pip-only users to import skasim, validate configuration, inspect CLI help, build docs, and run lightweight tests without Karabo installed. Full simulation execution should still require Karabo and fail with a clear runtime message when it is missing.

## Acceptance Criteria

- [ ] Importing lightweight skasim modules does not require Karabo.
- [ ] Running CLI help does not require Karabo.
- [ ] Configuration validation does not require Karabo.
- [ ] Simulation execution that needs Karabo fails with a clear installation message when Karabo is unavailable.
- [ ] Tests cover pip-only import/help behavior and missing-Karabo runtime failure.

## Blocked By

- Issue 01: Development Baseline

## Suggested Commit Message

`refactor: isolate Karabo behind runtime imports`
