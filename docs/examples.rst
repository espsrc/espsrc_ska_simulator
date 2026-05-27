Examples
========

Generated Sources
-----------------

Generate three point sources and create the default OSKAR dirty image::

    skasim --flux-density 1.0 5.0 10.0 \
      --telescope SKA1MID \
      --observation-time 600 \
      --frequency-mhz 1300 \
      --pixels 1024

Generated source polarization is optional. When supplied, the Stokes Q/U/V
lists must match the number of flux-density values::

    skasim --flux-density 1.0 5.0 \
      --stokes-q 0.1 0.2 \
      --stokes-u 0.0 0.1 \
      --stokes-v 0.0 -0.1 \
      --observation-time 600

Named-Catalog WSClean Smoke Check
-----------------------------------

This smoke-check command shape uses the MeerKAT telescope, the named MIGHTEE
catalog, the WSClean imager, and an explicit WSClean command::

    skasim --telescope MeerKAT \
      --catalog MIGHTEE \
      --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean" \
      --observation-time 60 \
      --frequency-mhz 1300 \
      --pixels 512 \
      --clean-iterations 100 \
      --output-dir smoke_mightee_wsclean

File-Backed FITS Smoke Check
----------------------------

This smoke-check command shape uses a FITS sky model file and the default
OSKAR dirty imager::

    skasim --model path/to/sources.fits \
      --telescope SKA1MID \
      --imager oskar-dirty \
      --observation-time 60 \
      --frequency-mhz 1300 \
      --pixels 512
