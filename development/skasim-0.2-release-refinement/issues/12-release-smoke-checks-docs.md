# Issue 12: Release Smoke Checks and Docs Polish

Type: AFK

## What to Build

Finalize the 0.2 release refinement by updating examples, validating the supported smoke command shape, and documenting the new CLI language. This slice should verify that the implementation remains aligned with the PRD.

## Acceptance Criteria

- [ ] README examples use named catalogue values and explicit imager values.
- [ ] User guide examples use source intensity language instead of unclear `--I` language.
- [ ] WSClean Singularity usage is documented through `--wsclean-command`.
- [ ] Smoke-check documentation includes a MeerKAT named-catalogue run and a file-backed FITS run shape.
- [ ] The local issue list and PRD remain consistent with implemented behavior.
- [ ] Tests and checks that can run in the available environment are documented.

## Blocked By

- Issue 03: Conda Environment and Install Docs
- Issue 11: Run Pipeline Records and Seconds Run IDs

## Suggested Commit Message

`docs: finalize skasim 0.2 release examples`
