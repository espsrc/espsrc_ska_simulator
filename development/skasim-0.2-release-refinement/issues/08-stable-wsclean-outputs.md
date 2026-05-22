# Issue 08: Stable WSClean Outputs

Type: AFK

## What to Build

Make WSClean output naming and collection run-scoped and stable. skasim should set a WSClean output prefix explicitly and collect only files belonging to the current image product.

## Acceptance Criteria

- [ ] WSClean runs use a stable run-specific `-name` prefix.
- [ ] WSClean output collection matches only the configured prefix.
- [ ] Output filenames are stable and do not encode scientific configuration parameters.
- [ ] Scientific and configuration metadata remain available in the manifest.
- [ ] Tests verify output discovery ignores unrelated WSClean files.

## Blocked By

- Issue 07: Configurable WSClean Command

## Suggested Commit Message

`refactor: make WSClean outputs run scoped`
