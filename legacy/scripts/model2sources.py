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
import numpy as np
from utils import show_exc



import xarray as xr
import os

if __name__ == "__main__":
    try:
        args = argparse.ArgumentParser()
        args.add_argument("--fits_model", type=str, required=True, help="Path to the FITS model file")
        args.add_argument("--output", type=str, help="Path to the output model file", default=None)
        args.add_argument("--sigma", type=float, help="Sigma threshold for source detection", default=10.0)
        args.add_argument("--freq", type=float, help="Frequency in MHz for source fluxes", default=1310)
        args.add_argument("--specindex", type=float, default=-0.7)
        args.add_argument("--components", type=int, help="Number of components to use for source modeling", default=None)
        args.add_argument("--rms", type=float, help="RMS noise level in Jy/beam to use for thresholding", default=None)
        args = args.parse_args()

        if os.path.exists(args.fits_model) is False:
            print (f"FITS model file {args.fits_model} does not exist!")
            sys.exit(1)



        f = args.fits_model
        img = SKAImage(f)
        if args.rms is not None:
            rms = args.rms
        else:
            rms = img.mad()
        if args.components is None:
            data = img.data2d
            threshold_value = args.sigma * rms
        else:
            sigma = img.mad()
            flat_data = img.data2d[~np.isnan(img.data2d)].flatten()    
            sorted_indices = np.argsort(flat_data)
            sorted_indices = sorted(sorted_indices, reverse=True)
            threshold_index = sorted_indices[args.components - 1]
            threshold_value = flat_data[threshold_index]
        # threshold_value /= (img.omega_beam() / img.omega_pix()) # Convert from Jy/beam to Jy/pixel

        img.tofits("temp_thresholded.fits")
        print(f"Using threshold value: {threshold_value:.6f}. Pixes above threshold: {np.sum(~np.isnan(img.data))}")
        oskar_sky = OskarSky.from_fits_file("temp_thresholded.fits", min_abs_val=threshold_value, frequency_hz = (args.freq * u.MHz).to(u.Hz).value, spectral_index = args.specindex)
        oskar_sky.filter_by_flux(min_flux_jy=threshold_value, max_flux_jy=np.inf)
        print (f"Total sources in OSKAR sky model: {oskar_sky.num_sources}")
        array_sky = oskar_sky.to_array()
        karabo_sky = SkyModel(sources=xr.DataArray(array_sky))

        print (f"Sources detected: {len(karabo_sky.sources)}")
        if (args.output is None):
            pickle_file = f.replace('.fits', '.kmod')
        else:
            if not args.output.endswith('.kmod'):
                pickle_file = f"{args.output}_kamod.kmod"
            else:
                pickle_file = args.output
        with open(pickle_file, 'wb') as pf:
            pickle.dump(karabo_sky, pf)
        print (f"Saved Karabo sky model to {pickle_file}")
    except Exception as e:
        print (f"Sorry, an error occurred.")
        print (show_exc(e))
        sys.exit(1)
