import os, sys
from datetime import datetime, timedelta
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.coordinates as acoord
from radio_beam import Beam
from astropy.table import Table
import xarray as xr

try:
    from karabo.simulation.telescope import Telescope # type: ignore
    from karabo.simulation.sky_model import SkyModel as KaraboSkyModel # type: ignore
except ImportError as e:
    print(f"Error importing Karabo modules: {e}")
import numpy as np
import json

def define_extra_units():
    # Definition of extra units
    u.def_unit('JY', 1*u.Jy)
    u.def_unit('DEG', 1*u.deg)
    u.def_unit('JY/BEAM', 1*u.Jy/u.sr)
    u.def_unit("HZ", 1*u.Hz)
    u.add_enabled_units(['JY', 'DEG', 'JY/BEAM', 'HZ'])

def mapping_unit(unit_str):

    unit_mapping = {
        'JY': 'Jy',
        'DEG': 'deg',
        'JY/BEAM': 'Jy/beam',
        'HZ': 'Hz'
    }
    if unit_str is None:
        return None
    return unit_mapping.get(unit_str, unit_str)


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)



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
                 major_axis = 0*u.arcsec, minor_axis = 0*u.arcsec, pa=0*u.deg,  true_redshift=0, obs_redshift=0, obj_id = None, resolved=False, isl_rms=0*u.Jy):
        # Initialize the source with its parameters, checking units

        list_of_units = [u.deg, u.deg, u.Jy, u.Jy, u.Jy, u.Jy, u.Hz, u.rad/(u.m**2), u.arcsec, u.arcsec, u.deg, u.Jy]
        list_of_values = [ra, dec, I, Q, U, V, ref_freq, rot_meas, major_axis, minor_axis, pa, isl_rms]
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
        self.resolved = resolved
        self.isl_rms = list_of_values[11]
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
            'resolved': self.resolved,
            'isl_rms': self.isl_rms.to(u.Jy).value
        }
    
    def from_array(array, colnames=['ra', 'dec', 'I', 'Q', 'U', 'V', 'ref_freq', 'spec_index', 'rot_meas', 'major_axis', 'minor_axis', 'pa', 'true_redshift', 'obs_redshift', 'resolved', 'isl_rms']):
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
        elif len(array) == 16:
            return Source(ra=array[ra_index] * u.deg, dec=array[dec_index] * u.deg, I=array[I_index] * u.Jy, 
                          Q=array[Q_index] * u.Jy, U=array[U_index] * u.Jy, V=array[V_index] * u.Jy,
                          ref_freq=array[ref_freq_index] * u.Hz, spec_index=array[spec_index_index], 
                          rot_meas=array[rot_meas_index] * (u.rad/(u.m**2)), 
                          major_axis=array[major_axis_index] * u.arcsec, minor_axis=array[minor_axis_index] * u.arcsec, 
                          pa=array[pa_index] * u.deg,
                          true_redshift=array[true_redshift_index],
                          obs_redshift=array[obs_redshift_index],
                          resolved=array[colnames.index('resolved')],
                          isl_rms=array[colnames.index('isl_rms')] * u.Jy)
        else:
            raise ValueError(f"Array must have 3, 6, or 12 elements (ra, dec, I, [Q, U, V], [ref_freq, spec_index, rot_meas, major_axis, minor_axis, pa, true_redshift, obs_redshift]). Your array has {len(array)} elements.")

    @staticmethod    
    def from_table_in_fits(table):
        # Generate sources from an astropy Table read from a FITS file
        col_equivs = [
            ['RA','ra'],
            ['DEC','dec'],
            ['STK_I','I', 'S_INT'],
            ['STK_Q','Q'],
            ['STK_U','U'],
            ['STK_V','V'],
            ['REFFREQ','ref_freq', 'NU_EFF'],
            ['SPECIDX','spec_index'],
            ['RM','rot_meas'],
            ['MAJ','major_axis', 'IM_MAJ'],
            ['MIN','minor_axis', 'IM_MIN'],
            ['PA','pa', 'IM_PA'],
            ['true_redshift','true_redshift'],
            ['obs_redshift','obs_redshift'],
            ['RESOLVED','resolved'],
            ['ISL_RMS','isl_rms']
        ]
        sources = []
        for alt_names in col_equivs:
            found = False
            for name in alt_names:
                if name in table.colnames:
                    colname = name
                    found = True
                    break
            if not found:
                colname = None
            alt_names.append(colname)  # The last element is the found column name or None
        for row in table:
            array = []
            for alt_names in col_equivs:
                colname = alt_names[-1]
                if colname is not None:
                    array.append(row[colname])
                else:
                    array.append(0)
            # array = [row['RA'], row['DEC'], row['STK_I']]
            # if 'STK_Q' in row.colnames:
            #     array.append(row['STK_Q'])
            # if 'STK_U' in row.colnames:
            #     array.append(row['STK_U'])
            # if 'STK_V' in row.colnames:
            #     array.append(row['STK_V'])
            # if 'REFFREQ' in row.colnames:
            #     array.append(row['REFFREQ'])
            # if 'SPECIDX' in row.colnames:
            #     array.append(row['SPECIDX'])
            # if 'RM' in row.colnames:
            #     array.append(row['RM'])
            # if 'MAJ' in row.colnames:
            #     array.append(row['MAJ'])
            # if 'MIN' in row.colnames:
            #     array.append(row['MIN'])
            # if 'PA' in row.colnames:
            #     array.append(row['PA'])
            # if 'true_redshift' in row.colnames:
            #     array.append(row['true_redshift'])
            # else:
            #     array.append(0)
            # if 'obs_redshift' in row.colnames:
            #     array.append(row['obs_redshift'])
            # else:
            #     array.append(0)
            # if 'RESOLVED' in row.colnames:
            #     array.append(row['RESOLVED'])
            # else:
            #     array.append(0)
            # if 'ISL_RMS' in row.colnames:
            #     array.append(row['ISL_RMS'])
            try:
                src = Source.from_array(array)
                if 'ID' in row.colnames:
                    src.obj_id = row['ID']
                sources.append(src)

            except Exception as e:
                print(show_exc(e))
                continue

        return sources

    def __str__(self):
        # Return a string representation of the source (only non-zero values)
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame='icrs')
        # Print RA/DEC in hms/dms format

        str2print = f"Source({coord.to_string('hmsdms')}, I={self.I:.6f})"

        json_values = self.to_json()
        for key, value in json_values.items():
            if key not in ['ra', 'dec', 'I'] and value != 0:
                value_with_unit = getattr(self, key)
                str2print += f", {key}={value_with_unit}"
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

    def get_flux(self, freq=None, alpha=None):
        # Return the total flux of the source
        if alpha is None:
            alpha = self.spec_index
        if freq is not None:
            # Calculate the flux at a given frequency using the spectral index
            flux = self.I * (freq / self.ref_freq) ** alpha
            return flux.to(u.Jy)
        return self.I.to(u.Jy)    
    
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
            
        def get_center(self, sources=None) -> SkyCoord:
            if self.phase_center is not None:
                return self.phase_center
            else:
                # Calculate the center of the sky model if phase_center is not set
                if sources is None:
                    sources = self.sources
                if sources.size > 0:
                    center_ra = np.mean(np.array(sources[:, 0])) * u.deg
                    center_dec = np.mean(np.array(sources[:, 1])) * u.deg
                    self.phase_center = SkyCoord(ra=center_ra, dec=center_dec, frame='icrs')
                    return self.phase_center
                else:
                    raise ValueError("SkyModel has no sources and phase_center is not set.")
                
        @staticmethod
        def from_fits_table(fits_file, log_file='sky_model.log', prefix='sky_model'):
            """        Load a SkyModel from a FITS table file. 
            Parameters:
            - fits_file: Path to the FITS file.
            Returns:
            - SkyModel object.
            """
            from astropy.io import fits
            if os.path.exists(fits_file):
                fits_table = Table.read(fits_file)
                sources = Source.from_table_in_fits(fits_table)
                sky_array = np.array([source.to_sky_model() for source in sources])
                skyModel = SkyModel(sky_array)
                skyModel.get_center()
                return skyModel
            else:
                raise FileNotFoundError(f"FITS file {fits_file} not found.")

except Exception as e:
    print(f"Error defining SkyModel class: {e}")