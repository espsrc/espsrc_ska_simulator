# Issue 11: Run Pipeline Records and Seconds Run IDs

Type: AFK

## What to Build

Deepen the run pipeline so a run consistently records milestones, structured outputs, errors, and weblog output while preserving the selected sky model source and image product behavior. Default run IDs should include seconds but no random suffix.

## Acceptance Criteria

- [ ] Default run IDs use `YYYYMMDD_HHMMSS_<telescope>`.
- [ ] No random suffix is added to default run IDs.
- [ ] The run pipeline builds the observation once for the run path rather than duplicating setup unnecessarily.
- [ ] Run records include structured outputs from sky model, simulation, imaging, manifest, log, and weblog phases where applicable.
- [ ] Failure paths save the manifest and render the weblog.
- [ ] Tests cover default run ID precision, success records, and failure records.

## Blocked By

- Issue 04: Named Sky Model Sources
- Issue 08: Stable WSClean Outputs
- Issue 10: Always-On Weblog

## Suggested Commit Message

`refactor: deepen run records and run ids`
