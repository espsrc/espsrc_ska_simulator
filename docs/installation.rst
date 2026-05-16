Installation
============

``skasim`` depends on the `Karabo <https://i4ds.github.io/Karabo-Pipeline/>`_
pipeline framework, the `OSKAR <https://ska-telescope.gitlab.io/sim/oskar/index.html>`_
simulation backend and, optionally, on the `WSClean <https://gitlab.com/aroffringa/wsclean>`_
imager for cleaned-image production.

Prerequisites
-------------

- Python ≥ 3.8
- `Karabo` with its full dependency stack (OSKAR, RASCIL, ska-sdp-datamodels)
- (Optional) ``wsclean`` binary on ``$PATH``, if you plan to use the cleaned imaging pathway

Start with a Miniconda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to use a clean Conda environment to manage the complex dependency stack::

    # Get Miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/miniconda3/bin/activate
    conda init bash
    conda install -n base conda-libmamba-solver

    # Create and activate environment
    conda create -n karabo python=3.9
    conda activate karabo

Install Karabo pipeline
~~~~~~~~~~~~~~~~~~~~~~~

Install the Karabo pipeline via Conda::

    conda install -c nvidia/label/cuda-11.7.0 -c i4ds -c conda-forge karabo-pipeline

Install RASCIL
~~~~~~~~~~~~~~

RASCIL must be installed from source. Ensure ``cmake`` and ``g++`` are available on your system::

    python -m pip install --upgrade pip setuptools wheel
    python -m pip install git+https://gitlab.com/ska-telescope/external/rascil.git

Install remaining dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install OSKAR and other required packages::

    python -m pip install --use-pep517 'git+https://github.com/OxfordSKA/OSKAR.git@master#egg=oskarpy&subdirectory=python'
    python -m pip install ska-sdp-datamodels --extra-index-url https://artefact.skao.int/repository/pypi-internal/simple
    python -m pip install thefuzz attrs radio-beam

Installing skasim
-----------------

Clone the repository and install in editable mode::

    git clone git@github.com:espsrc/espsrc_ska_simulator.git
    cd espsrc_ska_simulator
    pip install -e .

.. note::
   If you are installing into a pre-configured Karabo environment, use the ``--no-deps`` flag to avoid overwriting pinned dependencies:
   ``pip install -e . --no-deps``

WSClean (optional)
------------------

If you want to use the cleaned imaging pathway (``--cleaning`` flag), the
``wsclean`` binary must be installed separately and available on your
``$PATH``. On Debian/Ubuntu systems with the ``casacore`` stack::

    sudo apt install wsclean

Verifying the installation
--------------------------

After installation, confirm the CLI is available::

    skasim --help

You should see the help menu with options for telescope, frequency, imaging, and source settings.
