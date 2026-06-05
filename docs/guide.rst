User Guide
==========

This guide explains the main configuration language for a ``skasim`` run.

The Configuration Model
-----------------------

Every simulation is defined by a single ``SimConfig`` object, which nests:

- :class:`~skasim.config.ObsConfig` for observing setup;
- :class:`~skasim.config.ImgConfig` for image product settings;
- A list of :class:`~skasim.config.ModelEntry` objects for sky ingestion.

``ObsConfig`` validates the spectral grid. If bandwidth, channel count, and
channel width are all omitted, defaults are applied. If one value is omitted,
it is derived from the other two.

``ImgConfig`` records the selected imager. Multiple imaging passes are 
supported by supplying a list of ``ImgConfig`` entries, each with a unique 
``tag``::

    "imaging": [
        {"tag": "dirty", "imager": "oskar-dirty"},
        {"tag": "clean", "imager": "wsclean", "clean_iterations": 1000}
    ]

The ``--imager`` CLI flag overrides only the first entry. The default 
for a single pass is ``oskar-dirty``.

Sky Ingestion Models
--------------------

Version 0.2 introduced a typed models API. A run can include multiple model
entries of different types (though currently limited to one component-based
model per run).

Component Models
~~~~~~~~~~~~~~~~

**Type**: ``component_sky_model``

Used for component-based catalogs (Point sources, Gaussians).

Attributes:
- ``path``: Path to a FITS, JSON, pickle, or Karabo model file.
- ``catalog``: Named built-in catalog (e.g., ``MIGHTEE``, ``GLEAM``).
- ``flux_scale``: Global multiplier for source fluxes.

Example (JSON config)::

    "models": [
        {
            "type": "component_sky_model",
            "catalog": "MIGHTEE"
        }
    ]

Image Models
~~~~~~~~~~~~

Image-based models are injected via CASA ``ft`` into the ``MODEL_DATA`` column.

**1. Continuum I + Alpha** (``continuum_i_alpha``)

Injects a spatially extended source with a varying spectral index.

Attributes:
- ``stokes_i``: FITS image in Jy/pixel.
- ``alpha``: FITS image (spectral index map) on the same grid.
- ``reference_frequency_hz``: The frequency at which the ``stokes_i`` map is defined.

Semantics: :math:`I(\nu) = I_0 \cdot (\nu/\nu_0)^\alpha`.

**2. CASA Taylor Terms** (``casa_taylor_terms``)

Direct usage of existing CASA Taylor-term images (e.g., from a previous ``wsclean`` run).

Attributes:
- ``tt0``: Path to the CASA image table for Taylor term 0.
- ``tt1``: (Optional) Path to Taylor term 1.
- ``reference_frequency_hz``: The Taylor expansion reference frequency.

**3. Static Stokes Maps** (``static_stokes_maps``)

Injects static polarization maps.

Attributes:
- ``stokes_i``, ``stokes_q``, ``stokes_u``, ``stokes_v``: Paths to FITS images.

Spectral Reference Policy
-------------------------

When using image models, ``skasim`` automatically adjusts the model's spectral 
reference to the center of the observed band.

- **Idempotency**: The original model files are never modified. The pipeline 
  works on local copies in the ``work_dir``.
- **Correction**: If a spectral index is provided (via ``alpha`` map or ``tt1``), 
  pixel values are scaled using the power-law convention: 
  :math:`tt0' = tt0 \cdot (\nu_{new} / \nu_{old})^\alpha`.
- **Metadata**: The single-channel spectral coordinate (CRVAL4) of the local 
  CASA images is updated to the observation's center frequency.

Generated Source Fluxes (Legacy Mode)
--------------------------------------

If no typed models are provided, ``skasim`` can generate point sources using 
legacy CLI flags::

    skasim --flux-density 1.0 0.5 --observation-time 120

.. warning::
   Generated source flags are ignored if the ``models`` list is used.

WSClean Command
---------------

The WSClean command defaults to ``wsclean``. For containerized environments, 
pass the command explicitly::

    skasim --imager wsclean \
      --wsclean-command "singularity exec wsclean-3.10.sif wsclean"

Run Records
-----------

Each run writes a **Manifest** and an **Always-on Weblog**. 

.. note::
   See the :doc:`outputs` section for detailed information on run records.
