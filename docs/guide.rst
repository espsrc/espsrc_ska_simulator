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

Sky Model Sources
-----------------

Version 0.2 accepts one explicit sky model source per run:

1. a file-backed sky model supplied with ``--model``;
2. a named built-in catalogue supplied with ``--catalogue`` or ``--catalog``;
3. generated source intensities when no file or catalogue is supplied.

File-backed sky model
~~~~~~~~~~~~~~~~~~~~~

Supply a FITS, JSON, pickle, or Karabo model file::

    skasim --model my_sources.fits --center "05h00m00s -45d00m00s"

Named built-in catalogue
~~~~~~~~~~~~~~~~~~~~~~~~

Use catalogue names, not numeric IDs::

    skasim --catalogue MIGHTEE --center "05h00m00s -45d00m00s"
    skasim --catalog GLEAM --center "05h00m00s -45d00m00s"

Generated source intensities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If no file or catalogue is provided, ``skasim`` generates point sources around
the reference position. Each value in ``--flux-density`` creates one generated
source. ``--stokes-i`` is an alias with the same meaning::

    skasim --flux-density 1.0 0.5 0.2 --seconds 120
    skasim --stokes-i 1.0 0.5 0.2 --seconds 120

Source intensity flags are valid only in generated-source mode.

WSClean Command
---------------

The WSClean command defaults to ``wsclean``. On systems where WSClean is
provided through a wrapper or container, pass the command explicitly::

    skasim --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"

Run Records
-----------

Each run writes a manifest and an always-on weblog. The manifest contains
structured output records for logs, the manifest, visibility data, image
products, plots, and the weblog.

Output directories use the run id. By default, run ids use second precision::

    YYYYMMDD_HHMMSS_<telescope>
