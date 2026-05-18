# skasim — Spanish SRC SKA Simulator

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Package](https://img.shields.io/badge/pkg-0.1.0-orange)](https://github.com/espsrc/espsrc_ska_simulator)

**skasim** is a Python package for creating synthetic radio-interferometric observations of the SKA (Square Kilometre Array). It bridges user-supplied sky models and simulated imaging data products, supporting feasibility studies, definition of science cases and pipeline validation for the SKA community.

Born in the context of the CSIC4SKA project at IAA-CSIC, skasim wraps Karabo, OSKAR, and optionally WSClean behind a unified configuration model and a single CLI entrypoint.

```
Sky Model ──> OSKAR simulation ──> visibilities.MS ──> Image (OSKAR / WSClean)
```

---

## Key Features

- **Pydantic-validated configuration** — `SimConfig`, `ObsConfig`, `ImgConfig` catch parameter errors at construction time
- **Flexible sky-model inputs** — FITS catalogues, JSON source lists, Pickle/Karabo models, or built-in catalogues (MIGHTEE, GLEAM)
- **Two imaging pathways** — fast dirty imaging via OSKAR, or cleaned imaging via the external WSClean binary
---

## Quick Start

```bash

# Clone the repository
git clone https://github.com/espsrc/espsrc_ska_simulator.git
cd espsrc_ska_simulator

# Install (into a Karabo environment)
pip install -e .

# Run a simple simulation with 3 random point sources, SKA1-MID, clean imaging
skasim --I 1.0 5.0 10.0 --telescope SKA1MID --seconds 600 --freq 1300 --pixels 1024 --cleaning
```

See [full installation](https://github.com/espsrc/espsrc_ska_simulator/blob/main/docs/installation.rst) for environment setup.


---

## Dependencies

- **Python** ≥ 3.8
- **Karabo pipeline** (Conda-installable, includes OSKAR, RASCIL, ska-sdp-datamodels)
- **OSKAR** simulation backend (included in Karabo)
- **WSClean** (optional) — for cleaned imaging

---

## CLI Reference

```text
skasim --model <sky_model> --telescope <telescope> [options]

Key options:
  --freq MHz              Centre frequency
  --bandwidth MHz         Bandwidth
  --n-channels N          Number of channels
  --seconds N             Observation duration
  --cleaning              Use WSClean instead of OSKAR dirty
```

Run `skasim --help` for the complete reference, or check the full API documentation.

---

## Project Structure

```
src/skasim/
├── __init__.py         
├── cli.py              # argparse CLI
├── config.py           # Pydantic models: SimConfig, ObsConfig, ImgConfig
├── pipeline.py         # Orchestrator: sky → simulation → imaging
├── sky.py              # Source, SkyModel, catalogue loaders
├── imaging.py          # Dirty (OSKAR) / cleaned (WSClean) imaging
└── utils.py            # Constants, helpers, logger init
```

---

## References

1. **Karabo Pipeline** — I4DS framework for radio-astronomy simulations. [Documentation](https://i4ds.github.io/Karabo-Pipeline/) | [GitHub](https://github.com/i4ds/Karabo-Pipeline)
2. **OSKAR** — Oxford SKA Simulator. [Documentation](https://ska-telescope.gitlab.io/sim/oskar/) | [GitLab](https://gitlab.com/ska-telescope/sim/oskar)
3. **WSClean** — Widefield radio interferometric imager. [GitLab](https://gitlab.com/aroffringa/wsclean)
4. **RASCIL** — Radio Astronomy Simulation, Calibration and Imaging Library. [GitLab](https://gitlab.com/ska-telescope/external/rascil)
5. **ska-sdp-datamodels** — SKA Science Data Processor data models. [Repository](https://gitlab.com/ska-telescope/sdp/ska-sdp-datamodels)
