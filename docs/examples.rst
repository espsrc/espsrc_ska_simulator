Examples
========

Generated Sources
-----------------

Generate three point sources and create the default OSKAR dirty image::

    skasim --flux-density 1.0 5.0 10.0 \
      --telescope SKA1MID \
      --seconds 600 \
      --freq 1300 \
      --pixels 1024

The equivalent Stokes-I spelling is::

    skasim --stokes-i 1.0 5.0 10.0 --seconds 600

Named-Catalogue WSClean Smoke Check
-----------------------------------

This smoke-check command shape uses the MeerKAT telescope, the named MIGHTEE
catalogue, the WSClean imager, and an explicit WSClean command::

    skasim --telescope MeerKAT \
      --catalogue MIGHTEE \
      --imager wsclean \
      --wsclean-command "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean" \
      --seconds 60 \
      --freq 1300 \
      --pixels 512

File-Backed FITS Smoke Check
----------------------------

This smoke-check command shape uses a FITS sky model file and the default
OSKAR dirty imager::

    skasim --model path/to/sources.fits \
      --telescope SKA1MID \
      --imager oskar-dirty \
      --seconds 60 \
      --freq 1300 \
      --pixels 512
