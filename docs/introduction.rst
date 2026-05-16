Introduction
============

``skasim`` is a Python package for creating synthetic radio-interferometric
observations of the SKA (Square Kilometre Array). It bridges the gap between user-supplied
sky models and imaging data products.

Why a simulator?
----------------

The ``skasim`` simulator was born in the context of the CSIC4SKA project, aiming to support the CSIC community
in preparing for the SKAO Key Science Projects and devising competitive observational programs.

Simulated observations are essential for:

- **Feasibility studies** — can a given source be detected by SKA at the target frequency, bandwidth, and integration time?
- **Pipeline validation** — testing observation and imaging strategies before observing time is awarded.
- **Reproducibility** — a single ``skasim`` run captures the observational setup in a configuration file, making every
  simulation reproducible by design.

What skasim provides
--------------------

The goal of ``skasim`` is to deliver a low-friction wrapper that bundles existing astronomy tools to simulate observations and produce imaging products, facilitating the generation of synthetic images even for users with little experience on radio astronomy.

- **A unified configuration model** based on Pydantic — ``SimConfig``,
  ``ObsConfig``, ``ImgConfig`` — that validates parameters at construction time.
- **Multiple sky-model inputs**: FITS catalogues, JSON source lists, Pickle/Karabo models, and built-in catalogues (MIGHTEE, GLEAM).
- **Two imaging pathways**: fast dirty imaging via OSKAR, or cleaned (CLEAN-deconvolved) imaging via the external WSClean binary.
- **A single CLI entrypoint** — that exposes all options as command-line arguments.
- TODO: add more (web app, etc.)

How it works
-------------

::

    Sky Model ──> OSKAR simulation ──> visibilities.MS ──> Image
                    (Karabo)                              (OSKAR / WSClean)

The ``run()`` function orchestrates:

1. Telescope construction (loaded by Karabo)
2. Sky model loading or generation (from file, catalogue, or inline)
3. Observation definition (frequency, bandwidth, channels, duration)
4. Visibility simulation via the OSKAR backend
5. Imaging (dirty via OSKAR, or cleaned via WSClean)

Everything is controlled by a single ``SimConfig`` object, making the pipeline
reproducible and scriptable.
