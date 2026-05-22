# Issue 03: Conda Environment and Install Docs

Type: AFK

## What to Build

Add a supported conda environment for full simulation runs and update installation documentation to explain both pip-only and conda workflows. The conda environment should be named `skasim` and should start with Karabo as the full-runtime dependency path.

## Acceptance Criteria

- [ ] An `environment.yml` exists for the full simulation environment named `skasim`.
- [ ] Installation docs describe the pip-only workflow and what it supports.
- [ ] Installation docs describe the conda workflow and what it supports.
- [ ] Documentation does not promise a pip Karabo extra unless it has been verified.
- [ ] The WSClean Singularity command available on this machine is documented as an example, not as the default.

## Blocked By

- Issue 01: Development Baseline

## Suggested Commit Message

`docs: add pip and conda installation paths`
