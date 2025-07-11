import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.io import fits
from utils import show_exc


def guardar_distancia_fits(distance_matrix, wcs, filename):
    """
    Saves a distance matrix as a FITS file with WCS.

    Parameters:
    - distance_matrix: 2D array of distances in arcsec
    - wcs: corresponding WCS object
    - filename: output file name
    """
    hdu = fits.PrimaryHDU(data=distance_matrix, header=wcs.to_header())
    hdu.writeto(filename, overwrite=True)


def distance_matrix_in_arcsec(wcs, shape, ra_deg, dec_deg, radius_arcsec, pixscale=1*u.arcsec):
    """
    Returns a 2D matrix with angular distances in arcsec from the point (ra_deg, dec_deg).
    Pixels outside the given radius have value 0.

    Parameters:
    - wcs: astropy.wcs.WCS object
    - shape: (ny, nx) map size
    - ra_deg, dec_deg: center coordinates in degrees
    - radius_arcsec: maximum radius to consider, in arcsec

    Returns:
    - 2D distance matrix (in arcsec), with 0 outside the radius
    """
    try:
        ny, nx = shape
        if not isinstance(ra_deg, u.Quantity):
            ra_deg = ra_deg * u.deg

        if not isinstance(dec_deg, u.Quantity):
            dec_deg = dec_deg * u.deg

        if not isinstance(radius_arcsec, u.Quantity):
            radius_arcsec = radius_arcsec * u.arcsec

        if not isinstance(pixscale, u.Quantity):
            pixscale = pixscale * u.arcsec





        x_c, y_c = wcs.world_to_pixel(SkyCoord(ra=ra_deg, dec=dec_deg, frame='icrs'))
        # Apply to WCS
        wcs.wcs.cdelt = np.array([-1*pixscale.to(u.deg).value, pixscale.to(u.deg).value])  # degrees/pixel

        # Pixel scale (in degrees/pixel)
        no_dim = [i.to(u.deg).value for i in wcs.proj_plane_pixel_scales()]
        pixscale_deg = np.mean(np.abs(no_dim))
        radius_pix = radius_arcsec.to(u.deg).value / pixscale_deg  # Convert radius to pixels
        # Define window limits around the center
        x_min = int(max(0, np.floor(x_c - radius_pix)))
        x_max = int(min(nx, np.ceil(x_c + radius_pix)))
        y_min = int(max(0, np.floor(y_c - radius_pix)))
        y_max = int(min(ny, np.ceil(y_c + radius_pix)))

        # Create output matrix
        distance_matrix = np.zeros((ny, nx), dtype=np.float32)

        # Coordinates within the window
        y_indices, x_indices = np.mgrid[y_min:y_max, x_min:x_max]
        x_flat = x_indices.ravel()
        y_flat = y_indices.ravel()

        # Sky coordinates of those pixels
        world_coords = wcs.pixel_to_world(x_flat, y_flat)
        pix_coords = SkyCoord(world_coords.ra, world_coords.dec)

        # Center
        center = SkyCoord(ra=ra_deg, dec=dec_deg)

        # Distances
        separations = center.separation(pix_coords).arcsec * u.arcsec
        
        # Fill matrix where distance < radius
        mask = separations < radius_arcsec
        x_valid = x_flat[mask].astype(int)
        y_valid = y_flat[mask].astype(int)
        distances = separations[mask]

        distance_matrix[y_valid, x_valid] = distances
        return distance_matrix

    except Exception as e:
        print(show_exc(e))
        return np.zeros(shape, dtype=np.float32)
    
def setmask(matrix, r_in, r_out, value=0.0):
    """
    Masks values outside the range [r_in, r_out] and sets them to zero.
    
    Parameters:
    - matrix: 2D array of distances
    - r_in: inner radius in arcsec
    - r_out: outer radius in arcsec
    
    Returns:
    - Masked matrix with values outside the range set to zero.
    """
    mask = (matrix < r_in) | (matrix > r_out)
    matrix[mask] = value
    return matrix


if __name__ == "__main__":
    try:
        wcs = WCS(naxis=2)
        r_in = 0 * u.arcsec  # radio de 1 arcsec
        r_out = 30 * u.arcsec  # radio de 10 arcsec
        wcs.wcs.crpix = [50, 50]
        wcs.wcs.cdelt = np.array([-1/3600, 1/3600])
        wcs.wcs.crval = [150.0, 2.0]
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        shape = (100, 100)

        matriz = distance_matrix_in_arcsec(wcs, shape, 150.0, 2.0, r_out)  # radio 10 arcsec
        matriz = setmask(matriz, r_in.to(u.arcsec).value, r_out.to(u.arcsec).value, np.nan)



        mymask = np.isfinite(matriz)  # Mask for finite values
        mean_value = (np.nanmax(matriz) - np.nanmin(matriz)) / 2
        matriz[mymask] = -1 * np.abs(matriz[mymask] - mean_value)  # Normalizar a la mitad del rango
        matriz[mymask] += np.abs(np.min(matriz[mymask]))
        matriz[mymask] /= np.max(matriz[mymask])  # Normalizar entre 0 y 1

        matriz[(matriz < 0.95)] = 0

        # Extract the SkyCoords for pixes != 0
        non_zero_coords = np.argwhere(matriz != 0)
        sky_coords = wcs.pixel_to_world(non_zero_coords[:, 1], non_zero_coords[:, 0])
        print("SkyCoords of non-zero pixels:")
        for idx,coord in enumerate(sky_coords):
            print(idx, coord)





        guardar_distancia_fits(matriz, wcs, "angular_distance.fits")

        import matplotlib.pyplot as plt
        plt.imshow(matriz, origin='lower', cmap='viridis')
        plt.colorbar(label='Distancia angular (arcsec)')
        plt.title("Distancia angular desde el centro")
        plt.savefig("angular_distance.png")
    except Exception as e:
        print(show_exc(e))


