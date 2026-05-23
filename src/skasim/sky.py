import os
import sys
from datetime import datetime, timedelta

import astropy.coordinates as acoord
import numpy as np
import xarray as xr
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.table import Table
from astropy.time import Time
from astropy.utils.iers import conf as iers_conf
from loguru import logger

from .utils import show_exc
from radio_beam import Beam


class Source:
    def __init__(
        self,
        ra,
        dec,
        I,
        Q=0 * u.Jy,
        U=0 * u.Jy,
        V=0 * u.Jy,
        ref_freq=0 * u.Hz,
        spec_index=0,
        rot_meas=0 * u.rad / (u.m**2),
        major_axis=0 * u.arcsec,
        minor_axis=0 * u.arcsec,
        pa=0 * u.deg,
        true_redshift=0,
        obs_redshift=0,
        obj_id=None,
        resolved=False,
        isl_rms=0 * u.Jy,
    ):
        # initialize the source with its parameters, checking units

        list_of_units = [
            u.deg,
            u.deg,
            u.Jy,
            u.Jy,
            u.Jy,
            u.Jy,
            u.Hz,
            u.rad / (u.m**2),
            u.arcsec,
            u.arcsec,
            u.deg,
            u.Jy,
        ]
        list_of_values = [
            ra,
            dec,
            I,
            Q,
            U,
            V,
            ref_freq,
            rot_meas,
            major_axis,
            minor_axis,
            pa,
            isl_rms,
        ]
        for i, unit in enumerate(list_of_units):
            if not isinstance(list_of_values[i], u.Quantity):
                list_of_values[i] = list_of_values[i] * unit
            if list_of_values[i].unit != unit:
                try:
                    list_of_values[i] = list_of_values[i].to(unit)
                except u.UnitConversionError:
                    raise ValueError(
                        f"Value {list_of_values[i]} does not have the correct unit {unit}"
                    )
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
        self.obj_id = (
            obj_id if obj_id is not None else f"Source_{self.ra.value}_{self.dec.value}"
        )
        self.resolved = resolved
        self.isl_rms = list_of_values[11]
        self.coord = SkyCoord(
            ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame="icrs"
        )

    @staticmethod
    def from_name(name):
        try:
            source = acoord.get_icrs_coordinates(name)
        except NameResolveError:
            if name.upper() != "HCG16":
                raise
            source = SkyCoord(ra=32.390625 * u.deg, dec=-10.136389 * u.deg)
        if source is None:
            raise ValueError(f"Source {name} not found")
        return Source(source.ra, source.dec, 1 * u.Jy)

    def to_json(self, coords_fmt="deg"):
        # convert the source to a JSON serializable dictionary

        return {
            "ra": self.coords().to_string("hmsdms", sep=":").split()[0]
            if coords_fmt == "hmsdms"
            else self.coord.ra.to(u.deg).value,
            "dec": self.coords().to_string("hmsdms", sep=":").split()[1]
            if coords_fmt == "hmsdms"
            else self.coord.dec.to(u.deg).value,
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
            "resolved": self.resolved,
            "isl_rms": self.isl_rms.to(u.Jy).value,
        }

    def from_array(
        array,
        colnames=[
            "ra",
            "dec",
            "I",
            "Q",
            "U",
            "V",
            "ref_freq",
            "spec_index",
            "rot_meas",
            "major_axis",
            "minor_axis",
            "pa",
            "true_redshift",
            "obs_redshift",
            "resolved",
            "isl_rms",
        ],
    ):
        # create a Source object from a numpy array
        ra_index = colnames.index("ra")
        dec_index = colnames.index("dec")
        I_index = colnames.index("I")
        Q_index = colnames.index("Q") if "Q" in colnames else -1
        U_index = colnames.index("U") if "U" in colnames else -1
        V_index = colnames.index("V") if "V" in colnames else -1
        ref_freq_index = colnames.index("ref_freq") if "ref_freq" in colnames else -1
        spec_index_index = (
            colnames.index("spec_index") if "spec_index" in colnames else -1
        )
        rot_meas_index = colnames.index("rot_meas") if "rot_meas" in colnames else -1
        major_axis_index = (
            colnames.index("major_axis") if "major_axis" in colnames else -1
        )
        minor_axis_index = (
            colnames.index("minor_axis") if "minor_axis" in colnames else -1
        )
        pa_index = colnames.index("pa") if "pa" in colnames else -1
        true_redshift_index = (
            colnames.index("true_redshift") if "true_redshift" in colnames else -1
        )
        obs_redshift_index = (
            colnames.index("obs_redshift") if "obs_redshift" in colnames else -1
        )

        if len(array) < 3:
            raise ValueError("Array must have at least 3 elements (ra, dec, I)")
        if len(array) == 3:
            return Source(
                ra=array[ra_index] * u.deg,
                dec=array[dec_index] * u.deg,
                I=array[I_index] * u.Jy,
            )
        elif len(array) == 6:
            return Source(
                ra=array[ra_index] * u.deg,
                dec=array[dec_index] * u.deg,
                I=array[I_index] * u.Jy,
                Q=array[Q_index] * u.Jy,
                U=array[U_index] * u.Jy,
                V=array[V_index] * u.Jy,
            )
        elif len(array) == 12:
            return Source(
                ra=array[ra_index] * u.deg,
                dec=array[dec_index] * u.deg,
                I=array[I_index] * u.Jy,
                Q=array[Q_index] * u.Jy,
                U=array[U_index] * u.Jy,
                V=array[V_index] * u.Jy,
                ref_freq=array[ref_freq_index] * u.Hz,
                spec_index=array[spec_index_index],
                rot_meas=array[rot_meas_index] * (u.rad / (u.m**2)),
                major_axis=array[major_axis_index] * u.arcsec,
                minor_axis=array[minor_axis_index] * u.arcsec,
                pa=array[pa_index] * u.deg,
            )
        elif len(array) == 14:
            return Source(
                ra=array[ra_index] * u.deg,
                dec=array[dec_index] * u.deg,
                I=array[I_index] * u.Jy,
                Q=array[Q_index] * u.Jy,
                U=array[U_index] * u.Jy,
                V=array[V_index] * u.Jy,
                ref_freq=array[ref_freq_index] * u.Hz,
                spec_index=array[spec_index_index],
                rot_meas=array[rot_meas_index] * (u.rad / (u.m**2)),
                major_axis=array[major_axis_index] * u.arcsec,
                minor_axis=array[minor_axis_index] * u.arcsec,
                pa=array[pa_index] * u.deg,
                true_redshift=array[true_redshift_index],
                obs_redshift=array[obs_redshift_index],
            )
        elif len(array) == 16:
            return Source(
                ra=array[ra_index] * u.deg,
                dec=array[dec_index] * u.deg,
                I=array[I_index] * u.Jy,
                Q=array[Q_index] * u.Jy,
                U=array[U_index] * u.Jy,
                V=array[V_index] * u.Jy,
                ref_freq=array[ref_freq_index] * u.Hz,
                spec_index=array[spec_index_index],
                rot_meas=array[rot_meas_index] * (u.rad / (u.m**2)),
                major_axis=array[major_axis_index] * u.arcsec,
                minor_axis=array[minor_axis_index] * u.arcsec,
                pa=array[pa_index] * u.deg,
                true_redshift=array[true_redshift_index],
                obs_redshift=array[obs_redshift_index],
                resolved=array[colnames.index("resolved")],
                isl_rms=array[colnames.index("isl_rms")] * u.Jy,
            )
        else:
            raise ValueError(
                f"Array must have 3, 6, or 12 elements (ra, dec, I, [Q, U, V], [ref_freq, spec_index, rot_meas, major_axis, minor_axis, pa, true_redshift, obs_redshift]). Your array has {len(array)} elements."
            )

    @staticmethod
    def from_table_in_fits(table):
        # generate sources from an astropy Table read from a FITS file
        col_equivs = [
            ["RA", "ra"],
            ["DEC", "dec"],
            ["STK_I", "I", "S_INT"],
            ["STK_Q", "Q"],
            ["STK_U", "U"],
            ["STK_V", "V"],
            ["REFFREQ", "ref_freq", "NU_EFF"],
            ["SPECIDX", "spec_index"],
            ["RM", "rot_meas"],
            ["MAJ", "major_axis", "IM_MAJ"],
            ["MIN", "minor_axis", "IM_MIN"],
            ["PA", "pa", "IM_PA"],
            ["true_redshift", "true_redshift"],
            ["obs_redshift", "obs_redshift"],
            ["RESOLVED", "resolved"],
            ["ISL_RMS", "isl_rms"],
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
            alt_names.append(colname)
            # the last element is the found column name or None
        for row in table:
            array = []
            for alt_names in col_equivs:
                colname = alt_names[-1]
                if colname is not None:
                    array.append(row[colname])
                else:
                    array.append(0)

            try:
                src = Source.from_array(array)
                if "ID" in row.colnames:
                    src.obj_id = row["ID"]
                sources.append(src)

            except Exception as e:
                logger.error(show_exc(e))
                continue

        return sources


    def __str__(self):
        # return a string representation of the source (only non-zero values)
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.deg, u.deg), frame="icrs")
        # print RA/DEC in hms/dms format

        str2print = f"Source({coord.to_string('hmsdms')}, I={self.I:.6f})"

        json_values = self.to_json()
        for key, value in json_values.items():
            if key not in ["ra", "dec", "I"] and value != 0:
                value_with_unit = getattr(self, key)
                str2print += f", {key}={value_with_unit}"
        str2print += ")"
        return str2print

    def to_sky_model(self, reduced_form=False):
        # convert the source to a SkyModel object
        if reduced_form:
            return (self.ra.value, self.dec.value, self.I.value)
        else:
            return (
                self.ra.to(u.deg).value,
                self.dec.to(u.deg).value,
                self.I.to(u.Jy).value,
                self.Q.to(u.Jy).value,
                self.U.to(u.Jy).value,
                self.V.to(u.Jy).value,
                self.ref_freq.to(u.Hz).value,
                self.spec_index,
                self.rot_meas.value,
                self.major_axis.value,
                self.minor_axis.value,
                self.pa.value,
                self.true_redshift,
                self.obs_redshift,
            )

    def get_flux(self, freq=None, alpha=None):
        # return the total flux of the source
        if alpha is None:
            alpha = self.spec_index
        if freq is not None:
            # calculate the flux at a given frequency using the spectral index
            freq = u.Quantity(freq, u.Hz) if not isinstance(freq, u.Quantity) else freq
            flux = self.I * (freq / self.ref_freq) ** alpha
            return flux.to(u.Jy)
        return self.I.to(u.Jy)

    @property
    def flux(self):
        # return the total flux of the source
        return self.I

    def coords(self, frame="icrs"):
        # return the coordinates of the source
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
            date = datetime.now().strftime("%Y-%m-%d")

        coord = SkyCoord(ra=self.ra, dec=self.dec)
        location = EarthLocation(
            lat=telescope.centre_latitude * u.deg,
            lon=telescope.centre_longitude * u.deg,
            height=telescope.centre_altitude * u.m,
        )
        iers_conf.auto_download = False
        iers_conf.auto_max_age = None

        midnight = Time(f"{date} 00:00:00") + 12 * u.hour  # mediodía UTC
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

    @staticmethod
    def from_sky_model(data):
        """Reconstruct Source from 14-element tuple (inverse of to_sky_model)."""
        if len(data) == 3:
            return Source(ra=data[0], dec=data[1], I=data[2])
        return Source(
            ra=data[0], dec=data[1], I=data[2],
            Q=data[3], U=data[4], V=data[5],
            ref_freq=data[6], spec_index=data[7],
            rot_meas=data[8], major_axis=data[9],
            minor_axis=data[10], pa=data[11],
            true_redshift=data[12], obs_redshift=data[13],
        )

    def from_json(json_data):
        # create a Source object from a JSON dictionary
        return Source(
            ra=json_data["ra"] * u.deg,
            dec=json_data["dec"] * u.deg,
            I=json_data["I"] * u.Jy,
            Q=json_data["Q"] * u.Jy,
            U=json_data["U"] * u.Jy,
            V=json_data["V"] * u.Jy,
            ref_freq=json_data["ref_freq"] * u.Hz,
            spec_index=json_data["spec_index"],
            rot_meas=json_data["rot_meas"] * (u.rad / (u.m**2)),
            major_axis=json_data["major_axis"] * u.arcsec,
            minor_axis=json_data["minor_axis"] * u.arcsec,
            pa=json_data["pa"] * u.deg,
            true_redshift=json_data["true_redshift"],
            obs_redshift=json_data["obs_redshift"],
        )

    def to_fits_fmt(self):
        # convert the source to a FITS table format (dictionary)
        return {
            "RA": self.ra.to(u.deg).value,
            "DEC": self.dec.to(u.deg).value,
            "STK_I": self.I.to(u.Jy).value,
            "STK_Q": self.Q.to(u.Jy).value,
            "STK_U": self.U.to(u.Jy).value,
            "STK_V": self.V.to(u.Jy).value,
            "REFFREQ": self.ref_freq.to(u.Hz).value,
            "SPECIDX": self.spec_index,
            "RM": self.rot_meas.value,
            "MAJ": self.major_axis.to(u.arcsec).value,
            "MIN": self.minor_axis.to(u.arcsec).value,
            "PA": self.pa.to(u.deg).value,
            "true_redshift": self.true_redshift,
            "obs_redshift": self.obs_redshift,
            "RESOLVED": self.resolved,
            "ISL_RMS": self.isl_rms.to(u.Jy).value,
            "ID": self.obj_id,
        }


try:
    from karabo.simulation.sky_model import SkyModel as KaraboSkyModel
except ImportError:
    class KaraboSkyModel:
        """Small fallback for lightweight tests when Karabo is not installed."""

        def __init__(self, sources=None, **kwargs):
            self.sources = None if sources is None else np.asarray(sources)
            self.phase_center = None

        def add_point_sources(self, sources):
            sources_array = np.asarray(sources)
            if self.sources is None or len(self.sources) == 0:
                self.sources = sources_array
                return
            self.sources = np.vstack([self.sources, sources_array])


class SkyModel(KaraboSkyModel):
    phase_center = None

    def __init__(self, *args, **kwargs):
        # Call the parent constructor
        super().__init__(*args, **kwargs)
        if (self.sources is None) or (len(self.sources) == 0):
            return
        self.get_center()  # calculate the phase center if sources are provided

    def to_json(self):
        # convert the SkyModel to a JSON serializable list
        if self.sources is None:
            return []
        import xarray as xr
        if isinstance(self.sources, xr.DataArray):
            return [Source.from_sky_model(row.values).to_json() for row in self.sources]
        if isinstance(self.sources, np.ndarray):
            return [Source.from_sky_model(row).to_json() for row in self.sources]
        return [source.to_json() for source in self.sources]


    def show(self, **kwargs):
        if "block" not in kwargs:
            kwargs["block"] = False
        if "xlabel" not in kwargs:
            kwargs["xlabel"] = "RA (deg)"
        if "ylabel" not in kwargs:
            kwargs["ylabel"] = "DEC (deg)"
        logger.debug(self.phase_center)

        self.explore_sky(
            [
                self.phase_center.ra.to(u.deg).value,
                self.phase_center.dec.to(u.deg).value,
            ],
            **kwargs,
        )

    @staticmethod
    def from_json(json_data):
        # create a SkyModel object from a JSON list of sources
        try:
            sources = np.array(
                [Source.from_json(source).to_sky_model() for source in json_data]
            )
            skyModel = SkyModel(sources)
            center_ra = np.mean(sources[:, 0]) * u.deg
            center_dec = np.mean(sources[:, 1]) * u.deg

            skyModel.phase_center = SkyCoord(ra=center_ra, dec=center_dec, frame="icrs")
            return skyModel
        except Exception as e:
            print(e)
            return None

    @staticmethod
    def from_fits(
        fits_file,
        total_intensity=1 * u.Jy,
        fov=1 * u.deg,
        frequency=1 * u.GHz,
        log_file="sky_model.log",
        prefix="sky_model",
        t0=0,
    ):
        """Load a SkyModel from a FITS file.
        Parameters:
        - fits_file: Path to the FITS file.
        Returns:
        - SkyModel object.
        """
        import time

        import numpy as np
        from astropy.io import fits
        from astropy.wcs import WCS

        if os.path.exists(fits_file):
            source_ref = Source.from_name("HCG16")
            fits_data = fits.open(fits_file)
            fits_header = fits_data[0].header
            fits_data = fits_data[0].data
            img_pixels = int(fits_data.shape[2])
            fits_wcs = WCS(fits_header)
            sky_wcs = WCS(naxis=4)  # RA, DEC, Intensities, STOKES
            sky_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ", "STOKES"]
            sky_wcs.wcs.crpix = [fits_data.shape[2] // 2, fits_data.shape[3] // 2, 0, 0]
            sky_wcs.wcs.crval = [
                source_ref.ra.to(u.deg).value,
                source_ref.dec.to(u.deg).value,
                frequency.to(u.Hz).value,
                1.0,
            ]
            sky_wcs.wcs.cdelt = [
                fov.to(u.deg).value / img_pixels,
                fov.to(u.deg).value / img_pixels,
                1.0,
                1.0,
            ]
            sky_wcs.wcs.cunit = ["deg", "deg", "Hz", ""]

            fluxes = fits_data[0, 0, :, :]

            # check if total_intensity is a Quantity, if not, convert it
            if not isinstance(total_intensity, u.Quantity):
                total_intensity = total_intensity * u.Jy

            fluxes = (
                fluxes / np.max(fluxes) * total_intensity.to(u.Jy).value
            )  # Normalize to max intensity

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
            logger.debug("Starting conversion...")
            ra_list = []
            dec_list = []
            flux_list = []
            sum_weights = np.sum(fluxes)
            if total_intensity.to(u.Jy).value > 0:
                fluxes = fluxes / sum_weights * total_intensity.to(u.Jy).value
            else:
                logger.warning(
                    "Total_intensity is zero or negative, normalizing to 1 Jy"
                )
                fluxes = fluxes / sum_weights  # Normalize to 1 Jy

            for x, y in indices:
                world = sky_wcs.pixel_to_world(x, y, 0, 0)
                skycoord, freq, _ = world
                intensity = fluxes[x, y] * u.Jy
                ra_list.append(skycoord.ra.value)
                dec_list.append(skycoord.dec.value)
                flux_list.append(intensity.value)
                progress += 1
                logger.debug(
                    f"Progress: {progress:5.0f}/{total_pixels} ({progress / total_pixels * 100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds",
                    end="\r",
                )

                if progress in progress_to_print:
                    logger.debug(
                        f"Progress: {progress:5.0f}/{total_pixels} ({progress / total_pixels * 100:2.2f}%). Time elapsed: {time.time() - t0:.2f} seconds",
                    )

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
            # calculate the center of the sky model if phase_center is not set
            if sources is None:
                sources = self.sources
            if sources.size > 0:
                center_ra = np.mean(np.array(sources[:, 0])) * u.deg
                center_dec = np.mean(np.array(sources[:, 1])) * u.deg
                self.phase_center = SkyCoord(ra=center_ra, dec=center_dec, frame="icrs")
                return self.phase_center
            else:
                raise ValueError("SkyModel has no sources and phase_center is not set.")

    @staticmethod
    def from_fits_table(fits_file, log_file="sky_model.log", prefix="sky_model"):
        """Load a SkyModel from a FITS table file.
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
            return skyModel
        else:
            raise FileNotFoundError(f"FITS file {fits_file} not found.")
