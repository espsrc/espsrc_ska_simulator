import os, sys
from datetime import datetime, timedelta
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.coordinates as acoord
try:
    from karabo.simulation.telescope import Telescope # type: ignore
    from karabo.simulation.sky_model import SkyModel as KaraboSkyModel # type: ignore
except ImportError as e:
    print(f"Error importing Karabo modules: {e}")
import numpy as np



def printlog(fname, *args):
    # Print to console and log file
    print (f'[{datetime.now()}]',*args)
    with open(fname, 'a') as f:
        print(f'[{datetime.now()}]', *args, file=f)
        f.flush()
        os.fsync(f.fileno())
        f.close()

def show_exc(exception):
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    filename_rel = os.path.relpath(filename, os.path.dirname(__file__))
    app_folder = os.path.basename(os.path.dirname(__file__))
    return f'EXCEPTION IN ({filename_rel}:{lineno}): {exc_type} {exception} (APP: {app_folder})'


DIAMETERS = {
    'ALMA': 25 * u.m,
    'APEX': 12 * u.m,
    'ATCA': 22 * u.m,
    'CARMA': 10.4 * u.m,
    'GBT': 100 * u.m,
    'GMRT': 45 * u.m,
    'IRAM30M': 30 * u.m,
    'JCMT': 15 * u.m,
    'LOFAR': 25 * u.m,
    'MEERKAT': 13.5 * u.m,
    'MRT': 30 * u.m,
    'NRAO12M': 12 * u.m,
    'NRAO20M': 20 * u.m,
    'NRAO40M': 40 * u.m,
    'NRAO45M': 45 * u.m,
    'NRAO90M': 90 * u.m,
    'PARKES': 64 * u.m,
    'SMA': 6.5 * u.m,
    'SKA1LOW': 38 * u.m,
    'SKA1MID': 15 * u.m,
}

def get_diameter(telescope_name):
    """
    Returns the diameter of a telescope in meters.
    
    Parameters:
    - telescope_name: Name of the telescope (string).
    
    Returns:
    - Diameter in meters (astropy Quantity).
    """
    if telescope_name in DIAMETERS:
        return DIAMETERS[telescope_name]
    else:
        if ("SKA" in telescope_name or "SKA1" in telescope_name):
            # Handle SKA telescopes with specific names
            if "LOW" in telescope_name:
                return DIAMETERS['SKA1LOW']
            elif "MID" in telescope_name:
                return DIAMETERS['SKA1MID']
            else:
                raise ValueError(f"Telescope {telescope_name} not found. Available telescopes: {', '.join(DIAMETERS.keys())}")
        raise ValueError(f"Telescope {telescope_name} not found. Available telescopes: {', '.join(DIAMETERS.keys())}")

class Source:
    def __init__(self, ra, dec, I, Q=0 * u.Jy, U=0 * u.Jy, V=0 * u.Jy, ref_freq=0 * u.Hz, spec_index=0, rot_meas=0 * u.rad/(u.m**2), 
                 major_axis = 0*u.arcsec, minor_axis = 0*u.arcsec, pa=0*u.deg,  true_redshift=0, obs_redshift=0, obj_id = None):
        # Initialize the source with its parameters, checking units

        list_of_units = [u.deg, u.deg, u.Jy, u.Jy, u.Jy, u.Jy, u.Hz, u.rad/(u.m**2), u.arcsec, u.arcsec, u.deg]
        list_of_values = [ra, dec, I, Q, U, V, ref_freq, rot_meas, major_axis, minor_axis, pa]
        for i, unit in enumerate(list_of_units):
            if not isinstance(list_of_values[i], u.Quantity):
                list_of_values[i] = list_of_values[i] * unit
            if list_of_values[i].unit != unit:
                try:
                    list_of_values[i] = list_of_values[i].to(unit)
                except u.UnitConversionError:
                    raise ValueError(f"Value {list_of_values[i]} does not have the correct unit {unit}")
        self.ra = list_of_values[0]
        self.dec = list_of_values[1]
        self.I = list_of_values[2]
        self.Q = list_of_values[3]
        self.U = list_of_values[4]
        self.V = list_of_values[5]
        self.ref_freq = list_of_values[6]
        self.spec_index = spec_index
        self.rot_meas = list_of_values[7]
        self.major_axis = list_of_values[8]
        self.minor_axis = list_of_values[9]
        self.pa = list_of_values[10]
        self.true_redshift = true_redshift
        self.obs_redshift = obs_redshift
        self.obj_id = obj_id
        self.coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame='icrs')

    @staticmethod
    def from_name(name):
        source = acoord.get_icrs_coordinates(name)
        if source is None:
            raise ValueError(f"Source {name} not found")
        return Source(source.ra, source.dec, 1 * u.Jy)
    
    def to_json(self):
        # Convert the source to a JSON serializable dictionary
        return {
            "ra": self.ra.to(u.deg).value,
            "dec": self.dec.to(u.deg).value,
            "I": self.I.to(u.Jy).value,
            "Q": self.Q.to(u.Jy).value,
            "U": self.U.to(u.Jy).value,
            "V": self.V.to(u.Jy).value,
            "ref_freq": self.ref_freq.to(u.Hz).value,
            "spec_index": self.spec_index,
            "rot_meas": self.rot_meas.value,
            "major_axis": self.major_axis.to(u.arcsec).value,
            "minor_axis": self.minor_axis.to(u.arcsec).value,
            "pa": self.pa.to(u.deg).value,
            "true_redshift": self.true_redshift,
            "obs_redshift": self.obs_redshift,
        }
    
    def from_array(array, colnames=['ra', 'dec', 'I', 'Q', 'U', 'V', 'ref_freq', 'spec_index', 'rot_meas', 'major_axis', 'minor_axis', 'pa', 'true_redshift', 'obs_redshift']):
        # Create a Source object from a numpy array
        ra_index = colnames.index('ra')
        dec_index = colnames.index('dec')
        I_index = colnames.index('I')
        Q_index = colnames.index('Q') if 'Q' in colnames else -1
        U_index = colnames.index('U') if 'U' in colnames else -1
        V_index = colnames.index('V') if 'V' in colnames else -1
        ref_freq_index = colnames.index('ref_freq') if 'ref_freq' in colnames else -1
        spec_index_index = colnames.index('spec_index') if 'spec_index' in colnames else -1
        rot_meas_index = colnames.index('rot_meas') if 'rot_meas' in colnames else -1
        major_axis_index = colnames.index('major_axis') if 'major_axis' in colnames else -1
        minor_axis_index = colnames.index('minor_axis') if 'minor_axis' in colnames else -1
        pa_index = colnames.index('pa') if 'pa' in colnames else -1
        true_redshift_index = colnames.index('true_redshift') if 'true_redshift' in colnames else -1
        obs_redshift_index = colnames.index('obs_redshift') if 'obs_redshift' in colnames else -1

        if len(array) < 3:
            raise ValueError("Array must have at least 3 elements (ra, dec, I)")
        if len(array) == 3:
            return Source(ra=array[ra_index] * u.deg, dec=array[dec_index] * u.deg, I=array[I_index] * u.Jy)
        elif len(array) == 6:
            return Source(ra=array[ra_index] * u.deg, dec=array[dec_index] * u.deg, I=array[I_index] * u.Jy, 
                          Q=array[Q_index] * u.Jy, U=array[U_index] * u.Jy, V=array[V_index] * u.Jy)
        elif len(array) == 12:
            return Source(ra=array[ra_index] * u.deg, dec=array[dec_index] * u.deg, I=array[I_index] * u.Jy, 
                          Q=array[Q_index] * u.Jy, U=array[U_index] * u.Jy, V=array[V_index] * u.Jy,
                          ref_freq=array[ref_freq_index] * u.Hz, spec_index=array[spec_index_index], 
                          rot_meas=array[rot_meas_index] * (u.rad/(u.m**2)), 
                          major_axis=array[major_axis_index] * u.arcsec, minor_axis=array[minor_axis_index] * u.arcsec, 
                          pa=array[pa_index] * u.deg)
        elif len(array) == 14:
            return Source(ra=array[ra_index] * u.deg, dec=array[dec_index] * u.deg, I=array[I_index] * u.Jy, 
                          Q=array[Q_index] * u.Jy, U=array[U_index] * u.Jy, V=array[V_index] * u.Jy,
                          ref_freq=array[ref_freq_index] * u.Hz, spec_index=array[spec_index_index], 
                          rot_meas=array[rot_meas_index] * (u.rad/(u.m**2)), 
                          major_axis=array[major_axis_index] * u.arcsec, minor_axis=array[minor_axis_index] * u.arcsec, 
                          pa=array[pa_index] * u.deg,
                          true_redshift=array[true_redshift_index],
                          obs_redshift=array[obs_redshift_index])
        else:
            raise ValueError("Array must have 3, 6, or 12 elements (ra, dec, I, [Q, U, V], [ref_freq, spec_index, rot_meas, major_axis, minor_axis, pa, true_redshift, obs_redshift])")
    
    def __str__(self):
        # Return a string representation of the source (only non-zero values)
        str2print = f"Source(ra={self.ra}, dec={self.dec}, I={self.I}"

        json_values = self.to_json()
        for key, value in json_values.items():
            if key not in ['ra', 'dec', 'I'] and value != 0:
                str2print += f", {key}={value}"
        str2print += ")"
        return str2print
    
    def to_sky_model(self, reduced_form=False):
        # Convert the source to a SkyModel object
        if reduced_form:
            return (self.ra.value, self.dec.value, self.I.value)
        else:
            return( self.ra.to(u.deg).value, self.dec.to(u.deg).value, self.I.to(u.Jy).value, 
                self.Q.to(u.Jy).value, self.U.to(u.Jy).value, self.V.to(u.Jy).value,
                self.ref_freq.to(u.Hz).value, self.spec_index,
                self.rot_meas.value, self.major_axis.value,
                self.minor_axis.value, self.pa.value,
                self.true_redshift, self.obs_redshift)

    @property
    def flux(self):
        # Return the total flux of the source
        return self.I    
    
    def coords(self, frame='icrs'):
        # Return the coordinates of the source
        return SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame=frame)

    def get_best_observation_time(self, telescope, date=None):
        """
        Returns the local time at which an object with a given RA/Dec culminates (best observation time).
        
        Parameters:
        - ra_hours: Right Ascension in hours (float)
        - dec_degrees: Declination in degrees (float)
        - lat_deg: Observer's latitude in degrees (float)
        - lon_deg: Observer's longitude in degrees (float, positive to the East)
        - elevation_m: Altitude above sea level (optional)
        - date: Date as a string 'YYYY-MM-DD' (optional, defaults to today if not provided)
        - timezone_offset: Time difference relative to UTC (e.g., -6 for CDMX)

        Returns:
        - Best time.
        """

        if date is None:
            date = datetime.now().strftime('%Y-%m-%d')

        coord = SkyCoord(ra=self.ra, dec=self.dec)
        location = EarthLocation(lat=telescope.centre_latitude * u.deg, lon=telescope.centre_longitude * u.deg, height=telescope.centre_altitude * u.m)

        midnight = Time(f"{date} 00:00:00") + 12*u.hour  # Mediodía UTC
        delta = timedelta(minutes=1)

        best_time = None
        max_alt = -90

        for minutes in range(-360, 360):
            current_time = midnight + minutes * u.minute
            altaz = coord.transform_to(AltAz(obstime=current_time, location=location))
            if altaz.alt.deg > max_alt:
                max_alt = altaz.alt.deg
                best_time = current_time

        return best_time
    
    def from_json(json_data):
        # Create a Source object from a JSON dictionary
        return Source(
            ra=json_data['ra'] * u.deg,
            dec=json_data['dec'] * u.deg,
            I=json_data['I'] * u.Jy,
            Q=json_data['Q'] * u.Jy,
            U=json_data['U'] * u.Jy,
            V=json_data['V'] * u.Jy,
            ref_freq=json_data['ref_freq'] * u.Hz,
            spec_index=json_data['spec_index'],
            rot_meas=json_data['rot_meas'] * (u.rad/(u.m**2)),
            major_axis=json_data['major_axis'] * u.arcsec,
            minor_axis=json_data['minor_axis'] * u.arcsec,
            pa=json_data['pa'] * u.deg,
            true_redshift=json_data['true_redshift'],
            obs_redshift=json_data['obs_redshift']
        )

try:
    class SkyModel(KaraboSkyModel):
        phase_center = None

        def __init__(self, *args, **kwargs):
            # Call the parent constructor
            super().__init__(*args, **kwargs)
            if (self.sources is None) or (len(self.sources) == 0):
                return
            self.get_center()  # Calculate the phase center if sources are provided



        def to_json(self):
            # Convert the SkyModel to a JSON serializable dictionary
            return [source.to_json() for source in self.sources]
        
        def show(self, **kwargs):
            if "block" not in kwargs:
                kwargs["block"] = False
            if "xlabel" not in kwargs:
                kwargs["xlabel"] = "RA (deg)"
            if "ylabel" not in kwargs:
                kwargs["ylabel"] = "DEC (deg)"
            print(self.phase_center)
            
            self.explore_sky([self.phase_center.ra.to(u.deg).value, self.phase_center.dec.to(u.deg).value], **kwargs)
        
        @staticmethod
        def from_json(json_data):
            # Create a SkyModel object from a JSON list of sources
            try:
                sources = np.array([Source.from_json(source).to_sky_model() for source in json_data])
                skyModel = SkyModel(sources)
                center_ra = np.mean(sources[:, 0]) * u.deg
                center_dec = np.mean(sources[:, 1]) * u.deg

                skyModel.phase_center = SkyCoord(ra=center_ra, dec=center_dec, frame='icrs')
                return skyModel
            except Exception as e:
                print(show_exc(e))
                return None
        
        @staticmethod
        def from_fits(fits_file, total_intensity=1* u.Jy, fov=1 * u.deg, frequency=1 * u.GHz, log_file='sky_model.log', prefix='sky_model', t0=0):
            """        Load a SkyModel from a FITS file.
            Parameters:
            - fits_file: Path to the FITS file.
            Returns:
            - SkyModel object.
            """
            from astropy.io import fits
            from astropy.wcs import WCS
            import numpy as np
            import time

                
            if os.path.exists(fits_file):
                source_ref = Source.from_name("HCG16")
                fits_data = fits.open(fits_file)
                fits_header = fits_data[0].header
                fits_data = fits_data[0].data
                img_pixels = int(fits_data.shape[2])
                fits_wcs = WCS(fits_header)
                sky_wcs = WCS(naxis=4) # RA, DEC, Intensities, STOKES
                sky_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'FREQ', 'STOKES']
                sky_wcs.wcs.crpix = [fits_data.shape[2]//2, fits_data.shape[3] // 2, 0, 0]
                sky_wcs.wcs.crval = [source_ref.ra.to(u.deg).value, source_ref.dec.to(u.deg).value, frequency.to(u.Hz).value, 1.]
                sky_wcs.wcs.cdelt = [fov.to(u.deg).value / img_pixels, fov.to(u.deg).value / img_pixels, 1.0, 1.0]
                sky_wcs.wcs.cunit = ['deg', 'deg', 'Hz', '']

                fluxes = fits_data[0, 0, :, :] 

                # Check if total_intensity is a Quantity, if not, convert it
                if not isinstance(total_intensity, u.Quantity):
                    total_intensity = total_intensity * u.Jy

                fluxes = fluxes / np.max(fluxes) * total_intensity.to(u.Jy).value # Normalize to max intensity
                
                sources = []
                total_pixels = fluxes.size
                skyModel = SkyModel(wcs=sky_wcs)

                fluxes_nonzero = np.nonzero(fluxes)
                indices = np.array(fluxes_nonzero).T
                
                progress = 0
                total_pixels = indices.shape[0]
                progress_to_print = np.linspace(0, total_pixels, 11, dtype=int)
                max_flux = np.max(fluxes)
                t0 = time.time()
                printlog(log_file, "Starting conversion...")
                ra_list = []
                dec_list = []
                flux_list = []
                sum_weights = np.sum(fluxes)
                if total_intensity.to(u.Jy).value > 0:
                    fluxes = fluxes / sum_weights * total_intensity.to(u.Jy).value
                else:
                    printlog(log_file, "Warning: total_intensity is zero or negative, normalizing to 1 Jy")
                    fluxes = fluxes / sum_weights # Normalize to 1 Jy

                for x, y in indices:
                    world = sky_wcs.pixel_to_world(x, y, 0, 0)
                    skycoord, freq, _ = world 
                    intensity = fluxes[x, y] * u.Jy
                    ra_list.append(skycoord.ra.value)
                    dec_list.append(skycoord.dec.value)
                    flux_list.append(intensity.value)
                    progress += 1
                    print(f"Progress: {progress:5.0f}/{total_pixels} ({progress/total_pixels*100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds", end='\r')

                    if progress in progress_to_print:
                        printlog(log_file, f"Progress: {progress:5.0f}/{total_pixels} ({progress/total_pixels*100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds")
                
                np_samples = np.vstack((np.array(ra_list), np.array(dec_list))).transpose()
                np_fluxes = np.reshape(np.array(flux_list), (len(flux_list), 1))
                sky_array = np.hstack((np_samples, np_fluxes))
                skyModel = SkyModel(sky_array, wcs=sky_wcs)
                return skyModel
            else:
                raise FileNotFoundError(f"FITS file {fits_file} not found.")
            
        def get_center(self) -> SkyCoord:
            if self.phase_center is not None:
                return self.phase_center
            else:
                # Calculate the center of the sky model if phase_center is not set
                if self.sources.size > 0:
                    center_ra = np.mean(np.array(self.sources[:, 0])) * u.deg
                    center_dec = np.mean(np.array(self.sources[:, 1])) * u.deg
                    self.phase_center = SkyCoord(ra=center_ra, dec=center_dec, frame='icrs')
                    return self.phase_center
                else:
                    raise ValueError("SkyModel has no sources and phase_center is not set.")




                # pickle_path = os.path.join(os.path.dirname(__file__), f'{prefix}_sky_model.pkl')
                # with open(pickle_path, 'wb') as f:
                #     pickle.dump(skyModel.sources, f)
                # printlog (log_file, f"Sky model saved in {pickle_path}")

except Exception as e:
    print(f"Error defining SkyModel class: {e}")