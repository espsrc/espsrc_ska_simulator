# skasim 0.2 Release Refinement

This development folder contains the local planning artifacts for the skasim 0.2 release refinement.

- [PRD](./PRD.md)
- [Issue index and dependency graph](./ISSUES.md)
- [Local issues](./issues/)

## Baseline

- Current release baseline: `0.1.0` as declared in `pyproject.toml` and the public README badge.
- Target work: skasim `0.2` release refinement described by this folder.
- Domain language: the repository-level `CONTEXT.md` glossary is the source of truth for catalogue, sky model, sky model source, source intensity, run, image product, WSClean command, and weblog.

## Baseline Checks

- `python -m py_compile src/skasim/*.py` passes in the base environment.
- `python -m pytest -q` does not collect in the base environment because `pydantic` is missing and the local package is not installed on `PYTHONPATH`.

## Implementation Checks

- Lightweight test path used during implementation: `PYTHONPATH=/tmp/skasim-pydeps:src python -m pytest -q`.
- Syntax check used during implementation: `PYTHONPATH=/tmp/skasim-pydeps:src python -m py_compile src/skasim/*.py`.
- The temporary dependency path contains downloaded test dependencies for this session only; the supported full-runtime path remains `conda env create -f environment.yml`.

## Commit Message Rule

After implementing each issue, create a clear commit message that names the completed vertical slice. Each local issue includes a suggested commit message.
