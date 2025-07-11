import argparse
from astropy.io import fits
import numpy as np
import json
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import warnings
warnings.filterwarnings("ignore", category=UserWarning, append=True)
from scipy.optimize import curve_fit
import astropy.units as u

from utils import Source

def power_law(nu, S0, alpha, nu0):
    """
    Power law function for spectral index fitting.
    
    Parameters:
    - nu: frequency in Hz
    - S0: flux density at frequency nu0
    - alpha: spectral index
    - nu0: reference frequency in Hz
    """
    # Check if has units
    if isinstance(nu, u.Quantity):
        nu = nu.to(u.Hz).value
    if isinstance(nu0, u.Quantity):
        nu0 = nu0.to(u.Hz).value
    if isinstance(S0, u.Quantity):
        S0 = S0.to(u.Jy).value
    return S0 * (nu / nu0) ** (alpha)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spectral Analysis Script")
    parser.add_argument("cube", type=str, help="Path to the cube file for spectral analysis")
    parser.add_argument("json_srcs", type=str, help="Path to the JSON sources file")
    parser.add_argument("--index", type=int, default=0, help="Index of the source to process (default: 0)")

    
    args = parser.parse_args()

    cube_path = args.cube
    json_srcs_path = args.json_srcs

    cube = fits.open(cube_path)
    wcs = WCS(cube[0].header)


    json_data = json.load(open(json_srcs_path, 'r'))

    if args.index > 0:
        json_data = [json_data[args.index-1]]
    
    list_freqs = []
    list_spectra = []
    for json_source in json_data:
        source = Source.from_json(json_source)
        coords = SkyCoord(ra=source.ra, dec=source.dec, frame='icrs', unit='deg')
        x,y,_,_ = wcs.wcs_world2pix(coords.ra.to(u.deg).value, coords.dec.to(u.deg).value, (50*u.MHz).to(u.Hz).value,0,0)
        print(f"Processing source: {source} at pixel coordinates ({x}, {y}, {np.round(x)}, {np.round(y)})") 
        x, y = int(np.round(x)), int(np.round(y))

        # Extract the spectrum for the source (shape: (1, n_channels, 512, 512) : (Stoke, Channel, Y, X))
        frecuencies = cube[0].header['CRVAL3'] + np.arange(cube[0].header['NAXIS3']) * cube[0].header['CDELT3']
        frecuencies = frecuencies * u.Hz  # Convert to Hz
        spectrum = cube[0].data[0, :, y, x]  # Assuming the first Stokes parameter (I) and the given pixel
        list_freqs.append(frecuencies)
        list_spectra.append(spectrum)


    # Plot spectrum
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 5))
    # The x labels are in Hz. I want to show them in MHz
    for idx, (frecuencies, spectrum) in enumerate(zip(list_freqs, list_spectra)):
   
        plt.plot(frecuencies.to(u.MHz), spectrum, label=f'Spectrum at ({source.ra:.2f}, {source.dec:.2f})')
    plt.xlabel('Channel')
    plt.ylabel('Flux Density (Jy)')
    plt.title(fr'Spectrum from source $\alpha$={source.spec_index:.2f}')
    plt.legend()
    plt.savefig("spectrum_plot.png")


    try:
        spectrum[spectrum < 1e-6] = 1e-6  # Avoid negative values for logarithmic scale
        log_spectrum = np.log10(spectrum)  # Logarithmic scale for better visualization
        log_frecuencies = np.log10(frecuencies.to(u.Hz).value)  # Logarithmic scale for frequencies

        popt, pcov = curve_fit(power_law, frecuencies, spectrum)
        alpha, S0, nu0 = popt[1], popt[0], popt[2]
        print(f"Fitted spectral index: {alpha:.2f}, S0: {S0:.2f} Jy, nu0: {nu0:.2f} Hz")


    except Exception as e:
        print(f"Error fitting the spectral index: {e}")
        alpha, C = None, None
    # from scipy.stats import linregress
    # slope, intercept, r_value, p_value, std_err = linregress(log_frecuencies, log_spectrum)
    # print(f"Slope: {slope}, Intercept: {intercept}, R-squared: {r_value**2}, P-value: {p_value}, Std Err: {std_err}")

