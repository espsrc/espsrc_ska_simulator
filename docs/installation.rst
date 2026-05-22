Installation
============

``skasim`` supports two installation paths:

- a **pip-only** path for lightweight use, including imports, configuration
  validation, CLI help, documentation builds, and lightweight tests;
- a **conda full-runtime** path for executing simulations with Karabo, OSKAR,
  and the supported runtime stack.

Pip-only lightweight install
----------------------------

Use this path when you want to inspect the CLI, validate configuration, build
the documentation, or run tests that do not execute Karabo simulations::

    git clone git@github.com:espsrc/espsrc_ska_simulator.git
    cd espsrc_ska_simulator
    python -m pip install -e .

This install intentionally does not provide the full Karabo simulation runtime.
Commands such as ``skasim --help`` should work, but running a simulation requires
the conda environment below.

Conda full-runtime install
--------------------------

Use the checked-in environment file for simulation runs::

    conda env create -f environment.yml
    conda activate skasim

The environment is named ``skasim`` and installs ``karabo-pipeline`` from the
supported conda channels. It also installs this repository in editable mode.

WSClean
-------

The WSClean command defaults to ``wsclean``. If WSClean is available through a
wrapper or container, pass it explicitly with ``--wsclean-command``.

On this machine, the Singularity container can be used as an example::

    skasim --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"

This command is an environment-specific example, not the default.

Verifying The Installation
--------------------------

For the pip-only path, confirm that the CLI can be inspected::

    skasim --help

For the conda full-runtime path, activate ``skasim`` before running simulations.
