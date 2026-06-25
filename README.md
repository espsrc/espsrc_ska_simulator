<p align="center">
  <img src="docs/_static/skasim-logo.png" width="240" alt="skasim logo">
</p>

# skasim — Spanish SRC SKA Simulator

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Package](https://img.shields.io/badge/pkg-0.4.0-orange)](https://github.com/espsrc/espsrc_ska_simulator)

**skasim** is a Python package for creating synthetic radio-interferometric observations of the SKA (Square Kilometre Array). It bridges user-supplied sky models and simulated imaging data products, supporting feasibility studies, definition of science cases and pipeline validation for the SKA community.

Born in the context of the CSIC4SKA project at IAA-CSIC, skasim wraps Karabo, OSKAR, and optionally WSClean behind a unified configuration model and a single CLI entrypoint.

```
Sky Model ──> OSKAR simulation ──> visibilities.MS ──> Image (OSKAR / WSClean)
```

---

## Key Features

- **Pydantic-validated configuration** — `SimConfig`, `ObsConfig`, `ImgConfig` catch parameter errors and reject deprecated 0.1 fields
- **Typed sky-model API** — component catalogs, continuum I+α maps, CASA Taylor-term images, and static Stokes maps
- **Multi-imaging passes** — run several imager/tag combinations from a single visibility set
- **Two imaging pathways** — fast dirty imaging via OSKAR, or cleaned imaging via the external WSClean binary
- **Run records** — every run writes a structured manifest and an always-on weblog, including failed runs

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/espsrc/espsrc_ska_simulator.git
cd espsrc_ska_simulator

# Install (into a Karabo environment)
pip install -e .

# Run a small named-catalog simulation with WSClean imaging
skasim --output-dir smoke_mightee_wsclean \
  --telescope MeerKAT \
  --observation-time 60 \
  --frequency-mhz 1300 \
  --pixels 512 \
  --catalog MIGHTEE \
  --imager wsclean \
  --clean-iterations 100
```

See the [full installation guide](docs/installation.rst) and the [user guide](docs/guide.rst) for detailed environment setup and configuration options.

---

## Documentation

The full rendered documentation lives at **https://espsrc.github.io/espsrc_ska_simulator/**.
Realistic science demonstrations are collected in the public **science-case gallery** at https://espsrc.github.io/skasim-gallery/.

---

## Dependencies

- **Python** ≥ 3.10 (and < 3.13)
- **Karabo pipeline** (Conda-installable, includes OSKAR, RASCIL, ska-sdp-datamodels)
- **OSKAR** simulation backend (included in Karabo)
- **WSClean** (optional) — for cleaned imaging
- **python-casacore** / **CASA** (optional) — for Taylor-term and continuum image-model injection

---

## CLI Reference

```text
skasim --model <sky_model> --telescope <telescope> [options]

Key options:
  --config PATH               JSON configuration file
  --model PATH                Sky model file: FITS, JSON, pickle, or Karabo model
  --catalog NAME              Built-in catalog: MIGHTEE, GLEAM, or SKAMid
  --continuum-stokes-i PATH   Continuum image-model Stokes I FITS
  --continuum-alpha PATH      Continuum image-model spectral-index FITS
  --reference-frequency-hz N  Reference frequency for continuum image model
  --flux-density Jy...        Generated source Stokes I flux densities
  --stokes-q Jy...            Generated source Stokes Q values
  --stokes-u Jy...            Generated source Stokes U values
  --stokes-v Jy...            Generated source Stokes V values
  --telescope-version NAME    Array-assembly version (e.g. AA4, AA*)
  --show-telescopes           List available telescopes/versions and exit
  --frequency-mhz MHz         Central observing frequency
  --bandwidth-mhz MHz         Bandwidth
  --n-channels N              Number of channels
  --channel-width-mhz MHz     Channel width
  --observation-time SECONDS  Observation duration
  --center SEXAGESIMAL        Field centre, e.g. "10h01m35.1s 2d41m41s"
  --output-dir PATH           Exact output directory name
  --imager NAME               oskar-dirty (default) or wsclean
  --wsclean-command CMD       WSClean command or container invocation
  --clean-iterations N        WSClean CLEAN iterations
  --no-uv-coverage            Skip shadeMS UV-coverage plot
  --shadems-command CMD       shadeMS command/path
  --uv-coverage-canvas-size N Square canvas size for UV-coverage plot
  --rms                       Enable simple noise model
  --rms-value JY              RMS noise level (simple model)
  --noise-rms-start JY        OSKAR station noise RMS start
  --noise-rms-end JY          OSKAR station noise RMS end
  --overwrite                 Replace an existing output directory
```

Run `skasim --help` for the complete reference, or check the full API documentation.

Removed or renamed 0.1 options such as `--I`, `--stokes-i`,
`--fits`, `--json`, `--json-fg`, `--freq`, `--seconds`, `--prefix`, `--niter`,
`--scale-I`, `--cleaning`, numeric catalog IDs, and the Python config fields
`I`, `Q`, `U`, `V`, `cleaning`, `source_names`, `ref_freq_hz`, `json_fg`,
`output_prefix`, `niter`, `scale_I`, and `ImgConfig.algorithm` fail
intentionally in 0.2+.

---

## Project Structure

```
src/skasim/
├── __init__.py              # Public API: SimConfig, ObsConfig, ImgConfig
├── cli.py                   # argparse CLI
├── config.py                # Pydantic models: SimConfig, ObsConfig, ImgConfig, ModelEntry
├── pipeline.py              # Orchestrator: sky → simulation → imaging
├── sky.py                   # Source, SkyModel, catalog loaders
├── imaging.py               # Dirty (OSKAR) / cleaned (WSClean) imaging
├── loaders/                 # FITS catalog/image and CASA image-model loaders
├── manifest.py              # RunManifest, RunContext, milestone tracking
├── runtime.py               # Lazy Karabo/OSKAR/CASACore importers
├── weblog.py                # HTML report renderer
└── utils.py                 # Constants, helpers, logger init, shadeMS wrappers
```

---

## References

1. **Karabo Pipeline** — I4DS framework for radio-astronomy simulations. [Documentation](https://i4ds.github.io/Karabo-Pipeline/) | [GitHub](https://github.com/i4ds/Karabo-Pipeline)
2. **OSKAR** — Oxford SKA Simulator. [Documentation](https://ska-telescope.gitlab.io/sim/oskar/) | [GitLab](https://gitlab.com/ska-telescope/sim/oskar)
3. **WSClean** — Widefield radio interferometric imager. [GitLab](https://gitlab.com/aroffringa/wsclean)
4. **RASCIL** — Radio Astronomy Simulation, Calibration and Imaging Library. [GitLab](https://gitlab.com/ska-telescope/external/rascil)
5. **ska-sdp-datamodels** — SKA Science Data Processor data models. [Repository](https://gitlab.com/ska-telescope/sdp/ska-sdp-datamodels)
