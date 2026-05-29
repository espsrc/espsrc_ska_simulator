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

Use this procedure for simulation runs on the current development system. The
environment is named ``skasim`` and uses Python 3.10, Karabo, OSKAR, and WSClean
from conda packages. The editable pip install also installs Python dependencies
from ``pyproject.toml``, including ``shadems`` for UV-coverage plots::

    conda create -y -n skasim python=3.10
    conda install -y -n skasim --no-channel-priority \
      -c nvidia/label/cuda-11.7.0 \
      -c i4ds \
      -c conda-forge \
      karabo-pipeline "cuda-version=11.7"
    conda run -n skasim pip install -e .

The ``--no-channel-priority`` option is needed because strict channel priority
can block MPI/FFTW dependencies required by the Karabo runtime stack. The
``cuda-version=11.7`` pin is also needed: the OSKAR build installed with Karabo
expects CUDA 11 runtime libraries such as ``libcudart.so.11.0`` and
``libcufft.so.10``.

The checked-in ``environment.yml`` remains useful as a declaration of the
intended runtime dependencies, but the command sequence above is the verified
installation path on this machine.

WSClean
-------

The conda full-runtime install above installs WSClean through the
``karabo-pipeline`` dependency set. On the verified ``skasim`` environment,
``karabo-pipeline`` depends on ``wsclean`` and Conda installs
``wsclean 3.5.0`` from the ``i4ds`` channel.

WSClean checks whether OpenBLAS multi-threading is enabled and aborts if
``OPENBLAS_NUM_THREADS`` is not set to ``1``. The ``skasim`` pipeline sets this
environment variable automatically before launching WSClean, so pipeline runs
can succeed even when a direct shell command such as ``wsclean --version``
fails. To inspect WSClean directly, use::

    conda run -n skasim env OPENBLAS_NUM_THREADS=1 wsclean --version

The WSClean command defaults to ``wsclean``. If WSClean is available through a
wrapper or container instead, pass it explicitly with ``--wsclean-command``.

On this machine, the Singularity container can be used as an example::

    skasim --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"

This command is an environment-specific example, not the default.

shadeMS
-------

``skasim`` uses ``shadems`` to generate the UV-coverage plot displayed beside
the telescope layout in the weblog. ``shadems`` is installed from pip through
the package dependencies. To inspect it directly, use::

    conda run -n skasim shadems --help

The pipeline runs shadeMS with U and V axes, a square canvas, and writable
Matplotlib/Numba cache directories inside the run output directory.

Verifying The Installation
--------------------------

For the pip-only path, confirm that the CLI can be inspected::

    skasim --help

For the conda full-runtime path, verify the installed runtime before running
simulations::

    conda run -n skasim python -c "import karabo, oskar; print('karabo', karabo.__version__); print('oskar ok')"
    conda run -n skasim which oskar_sim_interferometer
    conda run -n skasim which wsclean
    conda run -n skasim which shadems
    conda run -n skasim skasim --help

``oskarpy`` is not required for the supported execution path; the successful
``import oskar`` check and a completed visibility simulation verify the OSKAR
backend used by Karabo.

A small WSClean smoke run with the named MIGHTEE catalog is::

    conda run -n skasim skasim \
      --output-dir smoke_mightee_wsclean_quick \
      --telescope MeerKAT \
      --observation-time 30 \
      --frequency-mhz 1300 \
      --bandwidth-mhz 25 \
      --n-channels 2 \
      --pixels 256 \
      --catalog MIGHTEE \
      --imager wsclean \
      --clean-iterations 20 \
      --overwrite

The verified smoke run completed with exit code 0 and wrote
``visibilities.MS``, a shadeMS UV-coverage plot, WSClean FITS products, PNG
previews, ``run_manifest.json``, and ``weblog.html`` under
``smoke_mightee_wsclean_quick/``.
