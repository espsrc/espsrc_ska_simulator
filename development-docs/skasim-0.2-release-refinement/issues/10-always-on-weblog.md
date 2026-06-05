# Issue 10: Always-On Weblog

Type: AFK

## What to Build

Ensure every run produces a weblog, including failed runs. The weblog is part of the run record and cannot be disabled.

## Acceptance Criteria

- [ ] Successful runs render a weblog.
- [ ] Failed runs render a weblog after the manifest is marked failed.
- [ ] The weblog displays run status, milestones, errors, and structured outputs.
- [ ] There is no CLI or configuration flag to disable weblog generation.
- [ ] Tests cover success and failure weblog rendering.

## Blocked By

- Issue 09: Structured Manifest Outputs

## Suggested Commit Message

`feat: render weblogs for every run`
