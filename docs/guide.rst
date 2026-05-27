User Guide
==========

This guide explains the main configuration language for a ``skasim`` run.

The Configuration Model
-----------------------

Every simulation is defined by a single ``SimConfig`` object, which nests:

- :class:`~skasim.config.ObsConfig` for observing setup;
- :class:`~skasim.config.ImgConfig` for image product settings.

``ObsConfig`` validates the spectral grid. If bandwidth, channel count, and
channel width are all omitted, defaults are applied. If one value is omitted,
it is derived from the other two.

``ImgConfig`` records the selected imager. The default is ``oskar-dirty``.
Use ``wsclean`` for cleaned imaging::

    skasim --imager wsclean --wsclean-command "wsclean"

Version 0.2 is a strict configuration boundary. Deprecated 0.1 Python fields
such as ``SimConfig.I``, ``SimConfig.Q``, ``SimConfig.U``, ``SimConfig.V``,
``SimConfig.cleaning``, ``SimConfig.source_names``, ``SimConfig.ref_freq_hz``,
``SimConfig.json_fg``, ``SimConfig.output_prefix``, ``SimConfig.niter``,
``SimConfig.scale_I``, and ``ImgConfig.algorithm`` are not accepted. Use
``source_flux_jy`` plus optional ``stokes_q_jy``, ``stokes_u_jy``, and
``stokes_v_jy`` for generated source fluxes and ``ImgConfig.imager`` for image
product selection.

Sky Model Sources
-----------------

Version 0.2 accepts one explicit sky model source per run:

1. a file-backed sky model supplied with ``--model``;
2. a named built-in catalog supplied with ``--catalog``;
3. generated source flux densities when no file or catalog is supplied.

File-backed sky model
~~~~~~~~~~~~~~~~~~~~~

Supply a FITS, JSON, pickle, or Karabo model file::

    skasim --model my_sources.fits --center "05h00m00s -45d00m00s"

Named built-in catalog
~~~~~~~~~~~~~~~~~~~~~~~~

Use catalog names, not numeric IDs::

    skasim --catalog MIGHTEE --center "05h00m00s -45d00m00s"

Generated source fluxes
~~~~~~~~~~~~~~~~~~~~~~~

If no file or catalog is provided, ``skasim`` generates point sources around
the reference position. Each value in ``--flux-density`` creates one generated
source::

    skasim --flux-density 1.0 0.5 0.2 --observation-time 120

Optional Stokes Q/U/V values can be supplied for generated sources. Each list
must match the ``--flux-density`` list length::

    skasim --flux-density 1.0 0.5 \
      --stokes-q 0.1 0.0 \
      --stokes-u 0.0 0.1 \
      --stokes-v 0.0 -0.1

Generated-source flux and polarization flags are valid only in generated-source mode.
If no generated-source intensity is supplied, the default generated source has
``10.0 Jy`` Stokes I flux density.

WSClean Command
---------------

The WSClean command defaults to ``wsclean``. On systems where WSClean is
provided through a wrapper or container, pass the command explicitly::

    skasim --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"

``--clean-iterations`` controls WSClean ``-niter``. The wrapper caps
WSClean ``-channels-out`` to the available simulated channel count, up to 8,
so small smoke runs can use fewer channels without invalid WSClean arguments.

Run Records
-----------

Each run writes a manifest and an always-on weblog. The manifest contains
structured output records for logs, the manifest, visibility data, image
products, plots, and the weblog.

Output directories use the run id. By default, run ids use second precision::

    YYYYMMDD_HHMMSS_<telescope>

Use ``--output-dir`` to choose the exact output directory name. Unlike the old
``--prefix`` behavior, the telescope name is not appended when ``--output-dir``
is supplied.

If an existing visibility output directory is present, ``skasim`` raises an
error unless ``--overwrite`` is supplied. There is no interactive overwrite
prompt in the run pipeline.
