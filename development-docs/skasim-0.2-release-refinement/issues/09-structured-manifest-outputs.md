# Issue 09: Structured Manifest Outputs

Type: AFK

## What to Build

Replace flat manifest output strings with structured output records that distinguish visibility data, image products, plots, logs, the manifest, and the weblog.

## Acceptance Criteria

- [ ] Manifest outputs are structured records with kind and path.
- [ ] Image product outputs include image product identity and imager metadata.
- [ ] Visibility, log, manifest, and weblog outputs can be identified by kind.
- [ ] Existing output information remains human-readable in serialized manifest JSON.
- [ ] Tests cover structured output serialization and backwards-readable manifest content.

## Blocked By

- Issue 02: Pip-Only Install and Lazy Karabo Runtime

## Suggested Commit Message

`feat: structure run manifest outputs`
