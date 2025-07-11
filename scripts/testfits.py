import oskar
import argparse
import numpy as np
from utils import Source
# from astropy.io import fits
# import astropy.units as u
from utils import show_exc
import json

def get_mad(data, iters=100, limit=0.99, stack=False, cutoff=2., mask=None):
    import warnings
    from astropy import stats
    try:
        warnings.simplefilter('ignore')
        try:
            if mask is not None:
                data[mask] = np.nan
        except Exception as e:
            show_exc ("Warning!!!! We can not masked the image")
        counter = 0
        next_mad = stats.mad_std(data, axis=None, ignore_nan=True)
        mad = 1e3
        mad_stack = []
        while (next_mad / mad < limit and counter < iters) or (counter == 0):
            mad = next_mad
            next_mad = stats.mad_std(data[np.where(data < cutoff * mad)], axis=None, ignore_nan=True)
            counter += 1
            mad_stack.append(mad)
            if next_mad <= 0:
                break

        if stack: 
            return(mad_stack)# * (u.Jy/beam)
        else:
            return(mad_stack[-1])# * (u.Jy/beam)
    except Exception as e:
        print(show_exc(e))


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Test FITS file reading for Sky model.")
        parser.add_argument("fits_file", type=str, help="Path to the FITS file to read")
        parser.add_argument("--sigma", type=float, default=3.0, help="Sigma threshold for source detection")
        args = parser.parse_args()
        sky = oskar.Sky.from_fits_file(args.fits_file)
        # Print some information about the sky model
        print(f"Sky model loaded from {args.fits_file}")
        sources = sky.to_array()
        mad = get_mad(sources[:, 2], iters=100, limit=0.99, stack=False, cutoff=2.0)
        sources = sources[(sources[:, 2] > args.sigma * mad)]

        json_data = []
        for source in sources:
            item = Source.from_array(source)
            json_data.append(item.to_json())

        json_file = args.fits_file.replace('.fits', '_sources.json')
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=4)
        print(f"Sources saved to {json_file}")
        
    




        # aux = SkyModel(sources)

        # sky_model = KaraboSkyModel(sources)
        # print(f"SkyModel created with {sky_model.num_sources} sources.")




        # sky_model = SkyModel(sources)
        # print(f"SkyModel created with {sky_model.num_sources} sources.")
        # sky_model.show(block=True)
    except Exception as e:
        print(show_exc(e))
