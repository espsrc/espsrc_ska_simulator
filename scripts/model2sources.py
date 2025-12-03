import oskar
import glob
import sys
from ska_img import SKAImage

from utils import Source, SkyModel
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import json
from oskar.sky import Sky as OskarSky
from karabo.simulation.sky_model  import SkyModel as KaraboSkyModel
import pickle
import argparse


# ref_freq = 1310 * u.MHz
# center = SkyCoord("10:00:27.4415357045 +2:20:57.0878906254", unit=(u.hourangle, u.deg), frame='icrs')
# model = SKAImage.new_model(center=center, fov=90*u.arcsec, pixels=1024, freq=ref_freq)

# fits_file = "./20251001_040_sources.fits"
# data = Table.read(fits_file)
# sources = Source.from_table_in_fits(data)
# print (f"Number of sources loaded: {len(sources)}")



# for src in sources:
#     coord_in_model = model.world_to_pixel(src.coord)
#     if coord_in_model is None:
#         print (f"Source {src.name} at {src.coord.to_string('hmsdms')} is out of image bounds.")
#         continue
#     (x,y) = coord_in_model
#     if x<0 or x>=model.size[0] or y<0 or y>=model.size[1]:
#         print (f"Source {src.name} at {src.coord.to_string('hmsdms')} is out of image bounds.")
#         continue
#     model.data[int(y), int(x)] += src.get_flux(ref_freq).to(u.Jy).value  # Add flux at 1 GHz

# model.tofits("test_model.fits")
# sys.exit()

import xarray as xr
import os

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("--fits_model", type=str, required=True, help="Path to the FITS model file")
    args.add_argument("--output", type=str, help="Path to the output model file", default=None)
    args.add_argument("--sigma", type=float, help="Sigma threshold for source detection", default=10.0)
    args.add_argument("--components", type=int, help="Number of components to use for source modeling", default=None)
    args = args.parse_args()

    if os.path.exists(args.fits_model) is False:
        print (f"FITS model file {args.fits_model} does not exist!")
        sys.exit(1)



    f = args.fits_model
    img = SKAImage(f)
    if args.components is None:
        img.data[(img.data < args.sigma * img.mad())] = 0.0
    else:
        # Keep only the N brightest components
        sigma = img.mad()
        flat_data = img.data.flatten()
        sorted_indices = flat_data.argsort()[::-1]  # Indices of sorted data in descending order
        threshold_index = sorted_indices[args.components - 1]
        threshold_value = flat_data[threshold_index]
        img.data[img.data < threshold_value] = 0.0
        print (f"Equivalent sigma: {threshold_value / sigma}")

    img.tofits("tmp.fits")
    oskar_sky = OskarSky.from_fits_file("tmp.fits")
    array_sky = oskar_sky.to_array()
    karabo_sky = SkyModel(sources=xr.DataArray(array_sky))
    print (f"Sources detected: {len(karabo_sky.sources)}")
    if (args.output is None):
        pickle_file = f.replace('.fits', '.karabo.mod')
    else:
        if not args.output.endswith('.karabo.mod'):
            pickle_file = f"{args.output}_kamod.karabo.mod"
        else:
            pickle_file = args.output
    with open(pickle_file, 'wb') as pf:
        pickle.dump(karabo_sky, pf)
    print (f"Saved Karabo sky model to {pickle_file}")
