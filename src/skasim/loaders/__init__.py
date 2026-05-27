"""skasim.loaders — Loaders that produce SkyModel from external sources.
"""

from .fits_catalogue import FitsCatalogLoader
from .fits_image import FitsImageLoader

__all__ = ["FitsCatalogLoader", "FitsImageLoader"]
