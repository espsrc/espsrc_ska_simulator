User Guide
==========

This guide explains the building blocks of a ``skasim`` simulation and how
to configure each one.

The configuration model
-----------------------

Every simulation is defined by a single ``SimConfig`` object, which nests
two additional sub-configurations:

- :class:`~skasim.config.ObsConfig` — observational parameters
  (frequency, bandwidth, channels, duration)
- :class:`~skasim.config.ImgConfig` — imaging parameters
  (pixels, field of view, Briggs robustness, algorithm)

These are Pydantic models with validation built in, so invalid values are
caught before the simulation starts.

ObsConfig — the observing setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from skasim.config import ObsConfig

    obs = ObsConfig(
        freq_mhz=1300.0,       # centre frequency (MHz)
        bandwidth_mhz=100.0,   # total bandwidth (MHz)
        n_channels=8,          # number of spectral channels
        seconds=600,           # integration time (s)
    )

Key points:

- ``freq_mhz`` defines the centre frequency (must be > 0).
- ``bandwidth_mhz``, ``n_channels`` and ``delta_freq_mhz`` form a
  *self-consistent triplet*. You must provide **at least two** of them;
  the third is derived automatically so that ``bandwidth = n_channels ×
  delta_freq``. If all three are omitted, the defaults ``100 MHz, 8
  channels, 12.5 MHz`` are used.
- If all three are provided, they must be mathematically consistent; otherwise a ``ValidationError`` is
  raised.
- ``seconds`` must be > 0.
- The phase centre is derived from the sky model automatically, or can be
  overridden via the ``center`` field in ``SimConfig``.

ImgConfig — imaging settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from skasim.config import ImgConfig

    img = ImgConfig(
        pixels=512,                 # image size (px, must be >= 64)
        fov_deg=None,               # field of view (deg); None = diffraction limit
        imaging_niter=1000,         # clean iterations
        robust=0.0,                 # Briggs robustness parameter
        algorithm="wsclean_clean",  # or "oskar_dirty"
    )

Key points:

- ``pixels`` must be ≥ 64.
- When ``fov_deg`` is ``None``, the field of view is computed from the
  telescope diameter and wavelength (diffraction limit, ∼1.2 λ/D).
- ``algorithm`` selects the imaging backend, but in the current version
  the CLI uses the ``--cleaning`` flag to switch; the pipeline reads
  ``SimConfig.cleaning`` rather than ``ImgConfig.algorithm``.

SimConfig — the main configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from skasim.config import SimConfig

    config = SimConfig(
        telescope="SKA1MID",
        observation=obs,
        imaging=img,
        # ... more fields below
    )

The top-level config holds:

**Telescope:**

- ``telescope`` — name string (e.g. ``"SKA1MID"``, ``"SKA1LOW"``,
  ``"MEERKAT"``).
- ``telescope_version`` — optional version string for telescopes that
  have multiple configurations.

**Sky model input (mutually exclusive resolution order):**

1. ``sky_file`` — path to a file (FITS, JSON, Pickle/Karabo model).
2. ``catalogue`` — integer selecting a built-in catalogue:
   ``1`` = MIGHTEE, ``2`` = GLEAM
3. Inline — if neither is given, a set of random point sources is
   generated around a reference position using the ``I`` list.

**Inline source parameters:**

- ``I`` — list of Stokes I intensities (Jy), one per source.
- ``Q``, ``U``, ``V`` — optional Stokes parameters.
- ``ref_freq_hz`` — reference frequency for spectral-index scaling.
- ``center`` — field centre as a string, e.g.
  ``"10h01m35.1s 2d41m41s"``.

**Noise:**

- ``rms`` — enable Karabo's built-in noise model (when ``True``, noise
  is modelled from the telescope receiver temperatures).

**Output:**

- ``output_prefix`` — directory name prefix for the run (default:
  timestamp-based).
- ``overwrite`` — whether to overwrite existing visibility files without
  prompting.

**Cleaning:**

- ``cleaning`` — if ``True``, uses WSClean instead of the OSKAR dirty
  imager.
- ``niter`` — number of CLEAN iterations (WSClean, default 5000).

Simulating sources
------------------

There are four ways to define the sky model. The pipeline tries them in
the order below.

1. From a file (recommended for repeatability)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Supply any file path to ``--model`` (or ``--fits`` / ``--json``)::

    skasim --model my_sources.json --center "05h00m00s -45d00m00s" --I 0.5

Supported formats:

- **FITS tables** — column names are matched automatically (RA, DEC,
  STK_I, STK_Q, STK_U, STK_V, REFFREQ, SPECIDX, etc.). Beam-unit
  conversion is attempted first; a pure-Jy fallback is used if beam
  units fail.
- **JSON** — a list of source dictionaries matching
  :meth:`Source.to_json <skasim.sky.Source.to_json>` output format.
  Supports all fields: ``ra``, ``dec``, ``I``, ``Q``, ``U``, ``V``,
  ``ref_freq``, ``spec_index``, ``rot_meas``, ``major_axis``,
  ``minor_axis``, ``pa``, ``true_redshift``, ``obs_redshift``,
  ``resolved``, ``isl_rms``.
- **Pickle / Karabo model** — a serialised ``SkyModel`` object
  (extensions ``.pkl``, ``.pickle``, ``.kmod``, ``.karabo.mod``).

2. From a built-in catalogue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``--catalogue N``::

    skasim --catalogue 1 --center "..." --I 0.5

Available catalogues (included with Karabo):

- ``1`` — MIGHTEE (radio continuum fields)
- ``2`` — GLEAM (Galactic and Extragalactic All-sky MWA)

3. Inline random sources
~~~~~~~~~~~~~~~~~~~~~~~~

If no file or catalogue is given, ``skasim`` generates random point
sources around a reference position (currently the galaxy cluster
HCG16). Each value in ``--I`` creates one source::

    skasim --I 1.0 0.5 0.2 --center "05h00m00s -45d00m00s" --seconds 120

4. From a Python script (advanced)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can construct the sky model programmatically::

    from skasim.sky import Source, SkyModel
    from astropy import units as u

    src = Source(
        ra=75.0 * u.deg,
        dec=-45.0 * u.deg,
        I=1.0 * u.Jy,
        spec_index=-0.7,
    )
    # Serialise to JSON for later re-use
    import json
    with open("my_source.json", "w") as f:
        json.dump([src.to_json()], f, indent=2)

The ``Source`` class handles unit conversion automatically — most
parameters accept either astropy ``Quantity`` objects or plain numbers
(with the appropriate default unit).

Running a simulation
--------------------

The simplest invocation::

    skasim --I 1.0 --seconds 60 \
           --ref-freq 700e6 --output-prefix my_test

This will:

1. Create a working directory ``my_test_<timestamp>_SKA1MID/``.
2. Generate one point source of 1 Jy at the phase center.
3. Simulate 60 seconds of SKA1-MID observations.
4. Produce a dirty image (PNG + FITS) inside the working directory.

Output structure
----------------

Each run creates a timestamped working directory containing:

::

    <prefix>_<telescope>/
    ├── <prefix>_<telescope>.log       # loguru log file
    ├── <prefix>_<telescope>_dirty.png # dirty / cleaned image (PNG)
    ├── <prefix>_<telescope>_dirty.fits
    ├── <prefix>_<telescope>_telescope.png  # telescope layout
    └── visibilities.MS/               # measurement set (directory)

When cleaning is enabled, additional image files are produced.

*work in progress*
