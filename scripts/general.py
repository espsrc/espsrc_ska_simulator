"""general.py — End-to-end interferometric simulation + imaging pipeline (Karabo/OSKAR/WSClean)

This script orchestrates a full “sky → visibilities → image” workflow aimed at SKA-style
interferometric simulations. It supports multiple sky-model ingestion paths (pickle, FITS
catalogue, JSON catalogue, built-in catalogues, or simple random sources), runs an OSKAR
visibility simulation through Karabo, and then produces either:

- a *dirty image* via Karabo's OSKAR imager; or
- a *cleaned image* via WSClean (invoked as an external binary through a custom command).

Outputs are written into a dedicated working directory named after `--prefix` (or a timestamp
default) plus the telescope identifier. The directory typically contains:

- `<prefix>.log` (plain-text log produced via `utils.printlog`)
- `<prefix>_visibilities.MS/` (Measurement Set with visibilities)
- `<prefix>_<telescope>_<version>_telescope.png` (telescope layout plot)
- Dirty imaging path:
  - `<prefix>_dirty.fits`, `<prefix>_dirty.png`
- Cleaning path (WSClean):
  - `wsclean-*.fits` intermediate images are renamed to include run parameters and
    corresponding `*.png` quicklooks are produced.

Assumptions / operational requirements
--------------------------------------
1) The Karabo framework is installed and configured, including OSKAR backends.
2) `wsclean` must be installed and discoverable on `$PATH` when `--cleaning` is used.
3) A local `utils.py` module must be available next to this script (or on `PYTHONPATH`),
   providing at least: `printlog`, `show_exc`, `Source`, `SkyModel`, `get_diameter`,
   `mapping_unit`, and unit helpers.

Astrophysical / interferometric conventions
-------------------------------------------
- The default field-of-view uses a diffraction-limited estimate $\mathrm{FoV} \approx 1.25\,\lambda/D$.
- Imaging cell size is set to `FoV / pixels` (in radians) and passed to imagers accordingly.
- Robust weighting and multiscale cleaning are controlled through the WSClean custom command.

Known sharp edges (documented but not changed)
----------------------------------------------
- `--I` is declared with `nargs='+'` but has a scalar default (`10`), which can lead to
  runtime errors when `--I` is omitted. Prefer `default=[10.0]`.
- FITS cube helper `create_fits_cube()` uses `np.concatenate(..., axis=1)` which likely does
  *not* create a standard spectral cube; `np.stack()` along a new axis is the more typical
  approach.

This file is intentionally kept close to the original runtime behavior; documentation and
comments were added without refactoring the control flow.

"""

try:
    from pandas.errors import AccessorRegistrationWarning
    import warnings
    warnings.filterwarnings(
        "ignore",
        category=AccessorRegistrationWarning,
)
except ImportError:
    pass

from matplotlib.colors import PowerNorm
import multiprocessing as mp
import json, sys, os, shutil, time, argparse, matplotlib, numpy as np
from datetime import timedelta, datetime
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
warnings.simplefilter('ignore', category=AstropyDeprecationWarning)
warnings.simplefilter("ignore", category=UserWarning)
import astropy.units as u
from typing import Dict
from astropy.units import UnitBase
from karabo.simulation.sky_model import SkyPrefixMapping, SkySourcesUnits
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
import karabo
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.imaging.imager_wsclean import (
    WscleanDirtyImager,
    WscleanImageCleaner,
    WscleanImageCleanerConfig,
    create_image_custom_command,
    TMP_PREFIX_CUSTOM,
    TMP_PURPOSE_CUSTOM,
    _WSCLEAN_BINARY,
    _get_command_prefix,
    )
from karabo.simulation.visibility import *
from karabo.simulator_backend import SimulatorBackend
from karabo.imaging.image import Image


from karabo.simulation.interferometer import InterferometerSimulation
from karabo.simulation.observation import Observation
# from karabo.simulation.sky_model import SkyModel
from karabo.simulation.telescope import Telescope
from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType

import argparse, pickle, re
from thefuzz import process, fuzz
from utils import printlog, show_exc, Source, get_diameter
import subprocess
import glob
from utils import SkyModel, define_extra_units, mapping_unit
# from karabo.simulation.sky_model import SkyModel



def letters_and_digits(s):
    """Return a string containing only letters and digits from the input string."""
    return ''.join(filter(str.isalnum, s)).lower().replace('versions', '')

def only_letters(s):
    """Return a string containing only letters from the input string."""
    return ''.join(filter(str.isalpha, s)).lower().replace('versions', '')

        

def iaa_create_image_custom_command(
    command: str,
    output_filenames: Union[str, List[str]] = "wsclean-image.fits",
) -> Union[Image, List[Image]]:
    """Create a dirty or cleaned image using your own command.

    Allows the use of the full WSClean functionality with all parameters.
    Command has to start with 'wsclean '.
    The working directory the command runs in will be a temporary directory.
    Use absolute paths to reference files or directories like the measurement set.

    Args:
        command: Command to execute.
        output_filenames: WSClean output filename(s)
            (relative to the working directory) that should be returned
            as Image objects. Can be a string for one file or a list of strings
            for multiple files.
            Example 1: "wsclean-image.fits"
            Example 2: ['wsclean-image.fits', 'wsclean-residual.fits']

    Returns:
        - If output_filenames is a **string**, returns an Image object of the file \
            output_filenames.
        - If output_filenames is a **list of strings**, returns a list of \
            Image objects, one object per filename in output_filenames.

    """

    try:
        tmp_dir = FileHandler().get_tmp_dir(
            prefix=TMP_PREFIX_CUSTOM,
            purpose=TMP_PURPOSE_CUSTOM,
        )
        expected_command_prefix = f"{_WSCLEAN_BINARY} "
        if not command.startswith(expected_command_prefix):
            raise ValueError(
                "Unexpected command. Expecting command to start with "
                f'"{expected_command_prefix}".'
            )
        # command = _get_command_prefix(tmp_dir) + command
        command = f"OPENBLAS_NUM_THREADS=1 {command}"
        print(f"WSClean command: [{command}]")
        completed_process = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            # Raises exception on return code != 0
            check=True,
        )
        print(f"WSClean output:\n[{completed_process.stdout}]")
        import glob



        tmp_files = glob.glob("wsclean-0*.fits")
        for tmp in tmp_files:
            try:
                os.remove(tmp)
            except Exception as e:
                print (show_exc(e))
        mfs_files = glob.glob("*-MFS-*.fits")
        print (f"Found MFS files: {mfs_files}")
        return [
            Image(path=mfs_file)
            for mfs_file in mfs_files
        ] if len(mfs_files) > 0 else None

    except Exception as e:
        printlog(f"{show_exc(e)}. \t\t Error running WSClean command: {command}")
        return None



def create_fits_cube(
    fits_files_pattern_or_list,
    output_filename="cube.fits",
    freq_header_key="CRVAL3",
    freq_unit="Hz"
):
    """
    Create a FITS cube by stacking 2D FITS images along a new frequency axis.

    Parameters:
        fits_files_pattern_or_list: Glob pattern (e.g., "*.fits") or list of FITS file paths.
        output_filename (str): Name of the output FITS cube.
        freq_header_key (str): Header keyword to extract frequency.
        freq_unit (str): Unit of the frequency axis.
    """
    import glob

    if isinstance(fits_files_pattern_or_list, str):
        fits_files = sorted(glob.glob(fits_files_pattern_or_list))
    else:
        fits_files = sorted(fits_files_pattern_or_list)

    if not fits_files:
        raise ValueError("No FITS files found.")

    data_list = []
    freqs = []

    for f in fits_files:
        with fits.open(f) as hdul:
            data_list.append(hdul[0].data)
            freq = hdul[0].header.get(freq_header_key)
            freqs.append(freq)

    cube_data = np.concatenate(data_list, axis=1)

    # Build header from first file
    with fits.open(fits_files[0]) as hdul:
        header = hdul[0].header.copy()

    

    header['NAXIS'] = 4
    header['NAXIS1'] = cube_data.shape[2]
    header['NAXIS2'] = cube_data.shape[1]
    header['NAXIS3'] = cube_data.shape[0]
    header['CTYPE3'] = 'FREQ'
    header['CUNIT3'] = freq_unit
    

    # Set CRVAL3 and CDELT3 if freqs are available
    if None not in freqs and len(set(freqs)) > 1:
        freqs = np.array(freqs)
        header['CRPIX3'] = 1
        header['CRVAL3'] = freqs[0]
        header['CDELT3'] = freqs[1] - freqs[0]  # Assumes linear spacing
    else:
        header['CRPIX3'] = 1
        header['CRVAL3'] = 1.0
        header['CDELT3'] = 1.0

    

    # Save cube
    hdu = fits.PrimaryHDU(data=cube_data, header=header)
    hdu.writeto(output_filename, overwrite=True)

    return output_filename

def get_telescope_version(telescope_name):
    """
    Get the version of the telescope from the name.
    If the name contains a version, return it, otherwise return None.
    """
    # Check if the telescope name contains a version
    version_telescope = None
    print (f"Getting version for telescope {telescope_name}")

    if telescope_name in get_args(OSKARTelescopesWithoutVersionType):
        print ("Telescope without version")
        return None


    if telescope_name in get_args(OSKARTelescopesWithVersionType):
        from karabo.simulation import telescope_versions
        available_versions = []
        available_objects = []
        for obj in dir(telescope_versions):
            # Check if the object is a class

            if isinstance(getattr(telescope_versions, obj), type):
                available_objects.append(getattr(telescope_versions, obj))
                available_versions.append(letters_and_digits(getattr(telescope_versions, obj).__name__))
        myversions = []
        similarity = process.extract(only_letters(telescope_name), available_versions, scorer=fuzz.ratio)
        for version, score in similarity:
            if score > 90:
                # get the index of the version in available_versions
                idx = available_versions.index(version)
                myversions.append(available_objects[idx])

        available_options = []
        available_parents = []
        for version_option in myversions:
            for option in version_option:
                available_parents.append(version_option)
                available_options.append(option)

        if len(available_options) == 0:
            print(f"No versions found for {telescope_name}. Exiting.")
            sys.exit(1)
        if len(available_options) == 1:
            version_telescope = available_options[0]
            print(f"Using version {version_telescope} for telescope {telescope_name}")
        elif len(myversions) > 1:
            
            option = None
            for idx, option in enumerate(available_options):
                print(f"{idx+1}. {option}")
            choice = input("Enter the number of the version you want to use: ")
            try:
                choice = int(choice)
                if choice < 1 or choice > len(myversions):
                    raise ValueError("Invalid choice")
                version_telescope = available_options[choice-1]
            except ValueError as e:
                print(f"Invalid choice: {e}. Exiting.")
                sys.exit(1)

        return version_telescope

def create_sky_model_from_file(fpath):
    # Examinate the file extension and then create the sky model
    ext = os.path.splitext(fpath)[-1].lower()
    if ext == '.pkl' or ext == '.pickle':
        with open(fpath, 'rb') as f:
            skyModel = pickle.load(f)
        return skyModel
    elif ext == '.fits' or ext == '.fit':
        printlog (log_file, f"Loading sources from FITS file {args.fits}")
        fits_path = fpath
        if not os.path.isabs(args.fits):
            fits_path = os.path.join(os.path.dirname(__file__), args.fits)

        if not os.path.isfile(fits_path):
            raise FileNotFoundError(f"File {fits_path} does not exist.")

        fits_file = fits.open(fits_path)

        unit_mapping: Dict[str, UnitBase] = {}


        cols_mapping = [int(i) for i in args.column_mapping.split(",")]
        str_cols = []
        for idx,col in enumerate(fits_file[1].columns):
            str_cols.append([idx,col.name])
            col_unit = mapping_unit(col.unit)
            unit_mapping[col.unit] = u.Unit(col_unit) if col_unit else u.dimensionless_unscaled
        prefix_mapping = SkyPrefixMapping(
            ra=fits_file[1].columns.names[cols_mapping[1]],
            dec=fits_file[1].columns.names[cols_mapping[2]],
            stokes_i=fits_file[1].columns.names[cols_mapping[3]],
            stokes_q=fits_file[1].columns.names[cols_mapping[4]] if cols_mapping[4] > -1 else None,
            stokes_u=fits_file[1].columns.names[cols_mapping[5]] if cols_mapping[5] > -1 else None,
            stokes_v=fits_file[1].columns.names[cols_mapping[6]] if cols_mapping[6] > -1 else None,
            spectral_index=fits_file[1].columns.names[cols_mapping[7]] if cols_mapping[7] > -1 else None,
            ref_freq=fits_file[1].columns.names[cols_mapping[8]] if cols_mapping[8] > -1 else None,
            rm=fits_file[1].columns.names[cols_mapping[9]] if cols_mapping[9] > -1 else None,
            major=fits_file[1].columns.names[cols_mapping[10]] if cols_mapping[10] > -1 else None,
            minor=fits_file[1].columns.names[cols_mapping[11]] if cols_mapping[11] > -1 else None,
            pa=fits_file[1].columns.names[cols_mapping[12]] if cols_mapping[12] > -1 else None,
            id=fits_file[1].columns.names[cols_mapping[0]] if cols_mapping[0] > -1 else None,
        )
        try:
            units_sources = SkySourcesUnits(
                stokes_i=u.Jy / u.beam,
                stokes_q=u.Jy / u.beam,
                stokes_u=u.Jy / u.beam,
                stokes_v=u.Jy / u.beam,
                ref_freq=u.MHz,
                major=u.arcsec,
                minor=u.arcsec,
                pa=u.deg,
                rm=u.rad / u.m**2,
            )
            skyModel = SkyModel.get_sky_model_from_fits(
                fits_file=fits_path,
                prefix_mapping=prefix_mapping,
                unit_mapping=unit_mapping,
                units_sources=units_sources,
                min_freq=None,
                max_freq=None,
                encoded_freq=None,
                memmap=False,
            )
        except u.core.UnitConversionError as e:
            printlog (log_file, f"Unit conversion error: {e}. Trying without beam units.")
            units_sources = SkySourcesUnits(
                stokes_i=u.Jy,
                stokes_q=u.Jy,
                stokes_u=u.Jy,
                stokes_v=u.Jy,
                ref_freq=u.MHz,
                major=u.arcsec,
                minor=u.arcsec,
                pa=u.deg,
                rm=u.rad / u.m**2,
            )
            skyModel = SkyModel.get_sky_model_from_fits(
                fits_file=fits_path,
                prefix_mapping=prefix_mapping,
                unit_mapping=unit_mapping,
                units_sources=units_sources,
                min_freq=None,
                max_freq=None,
                encoded_freq=None,
                memmap=False,
            )

        sources = skyModel.sources
        print(f"Loaded {len(sources)} sources from {fits_path}")
        return skyModel
    elif ext == '.json':
        printlog (log_file, f"Loading sources from JSON file {args.json}")
        if (not os.path.isabs(args.json)):
            fjson = os.path.join(os.path.dirname(__file__), args.json)
        else:
            fjson = args.json
        with open(fjson, 'r') as f:
            sources_data = json.load(f)
        sources = []
        for source_data in sources_data:
            source = Source.from_json(source_data)
            if args.scale_I is not None:
                source.I *= args.scale_I
            if source.ref_freq == 0:
                source.ref_freq = args.ref_freq[0] * u.Hz if args.ref_freq is not None else frequency.to(u.Hz)
            sources.append(source)
        if len(sources) == 0:
            raise ValueError("No sources found in JSON file")
        source_ref = sources[0]
        skyModel = SkyModel(sources=sources)
        return skyModel
    
 





if __name__ == "__main__":
    t0 = time.time()
    try:
        tmp_dir = FileHandler().get_tmp_dir(
            prefix=TMP_PREFIX_CUSTOM,
            purpose=TMP_PURPOSE_CUSTOM,
        )


        from karabo.simulation.telescope import OSKARTelescopesWithVersionType, OSKARTelescopesWithoutVersionType
        telescope_choices = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
        # Show the PID of the process
        print(f"PID: {os.getpid()}")

        # Use argparse for define two parameters; first parameter, a list with sources names; second parameters, show or write png
        argparser = argparse.ArgumentParser(description="SKAO simulation script")
        argparser.add_argument("--model", type=str, help="Sky model file", default=None)
        argparser.add_argument("--fits", type=str, help="FITS file with sources", default=None)
        argparser.add_argument("--json", type=str, help="JSON file with sources", default=None)
        argparser.add_argument("--json_fg", type=str, help="JSON file with foreground sources", default=None)
        argparser.add_argument("--center", type=str, help="Center of the field in RA,DEC (J2000), format: 10h01m35.1s 2d41m41s", default=None)
        argparser.add_argument("--prefix", type=str, help="Prefix for filenames", default=None)
        argparser.add_argument( "--telescope", choices=telescope_choices, type=str, help="Telescope to use for the simulation, default is 'SKA1LOW'", default="SKA-MID-AAstar", )
        argparser.add_argument("--I", type=float, default=10, help="Total Intensity in Jy (or max I if multiple sources)", nargs="+")
        argparser.add_argument("--Q", type=float, default=None, help="Q")
        argparser.add_argument("--U", type=float, default=None, help="U")
        argparser.add_argument("--V", type=float, default=None, help="V")
        argparser.add_argument("--ref_freq", type=float, nargs="+", default=None, help="Reference frequency in Hz")
        argparser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
        argparser.add_argument("--freq", type=float, help="Central Freq in MHz", default=200)
        argparser.add_argument("--bandwidth", type=float, help="Bandwidth in MHz", default=100)
        argparser.add_argument("--n_channels", type=int, help="Number of channels", default=0)
        argparser.add_argument("--delta_freq", type=float, help="Delta Freq in MHz", default=0)
        argparser.add_argument("--seconds", type=int, help="Observation Time in seconds", default=10)
        argparser.add_argument("--cleaning", action="store_true", help="Use cleaning algorithm")
        argparser.add_argument( "--pixels", type=int, help="Number of pixels in the image", default=512, )
        argparser.add_argument( "--imaging_niter", type=int, help="Number of iterations for the imager", default=1000, )
        argparser.add_argument("--fov", type=str, help="Field of view. Formats accepted: 1deg, 1arcsec, 1armin, 1rad. If not units, then assume degrees", default=None)
        argparser.add_argument("--robust", type=float, default=0.0, help="Robustness factor for imaging")
        argparser.add_argument("--show-telescopes", action="store_true", help="Show available telescopes and exit")

        argparser.add_argument("--rms", action="store_true", default=False, help="Enable RMS calculation")
        argparser.add_argument("--rms_value", type=float, default=0.0, help="RMS value in Jy")
        argparser.add_argument("--rms_sigma", type=float, default=3.0, help="Sigma factor for the RMS calculation to define the cleaning threshold")
        # Add noise_level argument
        # argparser.add_argument("--rms_start", type=float, default=0, help="Noise level start in Jy for RMS calculation")
        # # Add noise_level_var argument
        # argparser.add_argument("--rms_end", type=float, default=0, help="Noise level end in Jy for RMS calculation")
        argparser.add_argument("--niter", type=int, default=5000, help="Number of iterations for the imager")
        argparser.add_argument("--scale-I", type=float, default=1.0, help="Scale factor for I intensity")
        argparser.add_argument("--catalogue", type=int, default=0, help="Catalogue to use: 1 - MIGHTEE, 2 - GLEAM, 0 - JSON file or random sources")
        # Column mapping for FITS files format 0,1,2,3,4,5,6,7,8,9,10,11,12
        argparser.add_argument( "--column-mapping", type=str, help="Column mapping for FITS files format: id,ra,dec,stokes_i,stokes_q,stokes_u,stokes_v,spectral_index,ref_freq,rm,major,minor,pa,id", default="0,1,2,3,4,5,6,7,8,9,10,11,12", )
        args = argparser.parse_args()

        frequency = args.freq * u.MHz
        prefix = datetime.now().strftime("%Y%m%d_%H%M") if args.prefix is None else args.prefix

        if args.show_telescopes:
            print("Available telescopes:")
            for telescope in telescope_choices:
                print(f"- {telescope}")
            sys.exit(0)

        version_telescope = get_telescope_version(args.telescope)
        prefix = f"{prefix}_{args.telescope.replace('-','_')}" if args.telescope is not None else prefix
        os.makedirs(os.path.join(os.path.dirname(__file__), prefix), exist_ok=True)
        work_dir = os.path.join(os.path.dirname(__file__), prefix)
        os.chdir(work_dir)
        log_file = os.path.join(os.path.dirname(__file__), prefix, f'{prefix}.log')
        printlog (log_file, f"Command line: {' '.join(sys.argv)}")
        printlog (log_file, "System info")
        # Print command line
        printlog (log_file, f"\tKarabo version: {karabo.__version__}")
        try:
            total_memory_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
            total_memory_gb = total_memory_bytes / (1024. ** 3)
            printlog (log_file, f"\tTotal RAM: {total_memory_gb:.2f} GB")
        except Exception as e:
            printlog (log_file, f"\tError getting RAM memory (No LINUX OS?)")

        try:
            total_space_in_GB = shutil.disk_usage(os.path.expanduser(os.path.dirname(__file__))).total / (1024.0 ** 3)
            printlog (log_file, f"\tTotal disk space: {total_space_in_GB:.2f} GB")
        except Exception as e:
            printlog (log_file, f"\tError getting disk space)")

        printlog (log_file, f"\tNumber of CPU cores: {mp.cpu_count()}")

        backend = SimulatorBackend.OSKAR
        telescope_types = get_args(OSKARTelescopesWithVersionType) + get_args(OSKARTelescopesWithoutVersionType)
        if (args.telescope in telescope_types):
            if version_telescope is not None and version_telescope != args.telescope:
                printlog (log_file, f"Using telescope {args.telescope} with version {version_telescope}")
                telescope = Telescope.constructor(args.telescope, backend=backend, version=version_telescope)
            else:
                printlog (log_file, f"Using telescope {args.telescope} without version")
                telescope = Telescope.constructor(args.telescope, backend=backend)
                version_telescope = ''
        else:
            printlog (log_file, f"Telescope {args.telescope} not found. Loading SKA1LOW as default. Values accepted: {telescope_types}")
            telescope = Telescope.constructor("SKA1LOW", backend=backend)

        telescope.plot_telescope(file=f"{prefix}_{args.telescope}_{version_telescope}_telescope.png")    
        if args.fov is None:
            # Convert frequency to wavelength
            wavelength = frequency.to(u.m, equivalencies=u.spectral())
            fov = (1.25 * wavelength / get_diameter(args.telescope.upper())) * u.rad
        else:
            try:
                fov = u.Quantity(args.fov)

                if fov.unit.is_equivalent(u.deg):
                    fov = fov.to(u.rad)
                elif fov.unit.is_equivalent(u.rad):
                    fov = fov.t(u.rad)
                else:
                    fov = fov.value * u.deg
            except Exception as e:
                print(show_exc(e))
                fov = float(args.fov) * u.deg
                fov = fov.to(u.rad)

        skyModel = None
        skyModel_fg = None



        try:
            if args.model is not None:
                # printlog (log_file, f"Loading sky model from file {args.model}")
                # model_path = args.model
                # if not os.path.isabs(args.model):
                #     model_path = os.path.join(os.path.dirname(__file__), args.model)
                # if not os.path.isfile(model_path):
                #     raise FileNotFoundError(f"Sky model file {model_path} does not exist.")
                # with open(model_path, 'rb') as mf:
                #     skyModel = pickle.load(mf)
                # sources = skyModel.sources
                # printlog (log_file, f"Loaded {len(sources)} sources from sky model {model_path}")
                skyModel = create_sky_model_from_file(args.model)
                sources = skyModel.sources
                printlog (log_file, f"Loaded {len(sources)} sources from sky model {args.model}")
            else:
                if args.catalogue > 0:
                    if args.catalogue == 1:
                        printlog (log_file, "Loading sources from MIGHTEE catalogue")
                        skyModel = SkyModel.get_MIGHTEE_Sky()
                    elif args.catalogue == 2:
                        printlog (log_file, "Loading sources from GLEAM catalogue")
                        skyModel = SkyModel.get_GLEAM_Sky()
                    elif args.catalogue == 3:
                        # Check if exists ./SKAMid_B1_8h_v3.fits; else download from https://owncloud.ia2.inaf.it/index.php/s/H1rAR0A9qmXBbB5/download
                        skamid_path = os.path.join(os.path.dirname(__file__), 'SKAMid_B1_8h_v3.fits')
                        if os.path.exists(skamid_path):
                            printlog (log_file, f"Loading sources from SKAMid catalogue {skamid_path}")
                            skyModel = SkyModel.get_sky_model_from_fits(fits_file=skamid_path)
                        else:
                            printlog (log_file, f"Catalogue {skamid_path} not found. Downloading from https://owncloud.ia2.inaf.it/index.php/s/H1rAR0A9qmXBbB5/download")
                            skamid_url = "https://owncloud.ia2.inaf.it/index.php/s/H1rAR0A9qmXBbB5/download"
                            import requests
                            response = requests.get(skamid_url)
                            if response.status_code == 200:
                                with open(skamid_path, 'wb') as f:
                                    f.write(response.content)
                                printlog (log_file, f"Catalogue {skamid_path} downloaded successfully")
                                skyModel = SkyModel.get_sky_model_from_fits(fits_file=skamid_path)
                            else:
                                printlog (log_file, f"Error downloading catalogue {skamid_url}. Status code: {response.status_code}")
                                sys.exit(1)
                    else:
                        printlog (log_file, f"Catalogue {args.catalogue} not found. Available catalogues: 1 - MIGHTEE, 2 - GLEAM")
                        sys.exit(1)

                    sources = skyModel.sources
                elif args.fits is not None:
                    # printlog (log_file, f"Loading sources from FITS file {args.fits}")
                    # fits_path = args.fits
                    # if not os.path.isabs(args.fits):
                    #     fits_path = os.path.join(os.path.dirname(__file__), args.fits)

                    # if not os.path.isfile(fits_path):
                    #     raise FileNotFoundError(f"File {fits_path} does not exist.")
                    # # survey_file = fits.open(args.fits)

                    # fits_file = fits.open(fits_path)

                    # unit_mapping: Dict[str, UnitBase] = {}

                    # cols_mapping = [int(i) for i in args.column_mapping.split(",")]
                    # str_cols = []
                    # for idx,col in enumerate(fits_file[1].columns):
                    #     str_cols.append([idx,col.name])
                    #     col_unit = mapping_unit(col.unit)
                    #     unit_mapping[col.unit] = u.Unit(col_unit) if col_unit else u.dimensionless_unscaled
                    # prefix_mapping = SkyPrefixMapping(
                    #     ra=fits_file[1].columns.names[cols_mapping[1]],
                    #     dec=fits_file[1].columns.names[cols_mapping[2]],
                    #     stokes_i=fits_file[1].columns.names[cols_mapping[3]],
                    #     stokes_q=fits_file[1].columns.names[cols_mapping[4]] if cols_mapping[4] > -1 else None,
                    #     stokes_u=fits_file[1].columns.names[cols_mapping[5]] if cols_mapping[5] > -1 else None,
                    #     stokes_v=fits_file[1].columns.names[cols_mapping[6]] if cols_mapping[6] > -1 else None,
                    #     spectral_index=fits_file[1].columns.names[cols_mapping[7]] if cols_mapping[7] > -1 else None,
                    #     ref_freq=fits_file[1].columns.names[cols_mapping[8]] if cols_mapping[8] > -1 else None,
                    #     rm=fits_file[1].columns.names[cols_mapping[9]] if cols_mapping[9] > -1 else None,
                    #     major=fits_file[1].columns.names[cols_mapping[10]] if cols_mapping[10] > -1 else None,
                    #     minor=fits_file[1].columns.names[cols_mapping[11]] if cols_mapping[11] > -1 else None,
                    #     pa=fits_file[1].columns.names[cols_mapping[12]] if cols_mapping[12] > -1 else None,
                    #     id=fits_file[1].columns.names[cols_mapping[0]] if cols_mapping[0] > -1 else None,
                    # )
                    # units_sources = SkySourcesUnits(
                    #     stokes_i=u.Jy / u.beam,
                    #     stokes_q=u.Jy / u.beam,
                    #     stokes_u=u.Jy / u.beam,
                    #     stokes_v=u.Jy / u.beam,
                    #     ref_freq=u.MHz,
                    #     major=u.arcsec,
                    #     minor=u.arcsec,
                    #     pa=u.deg,
                    #     rm=u.rad / u.m**2,
                    # )
                    # skyModel = SkyModel.get_sky_model_from_fits(
                    #     fits_file=fits_path,
                    #     prefix_mapping=prefix_mapping,
                    #     unit_mapping=unit_mapping,
                    #     units_sources=units_sources,
                    #     min_freq=None,
                    #     max_freq=None,
                    #     encoded_freq=None,
                    #     memmap=False,
                    # )
                    # # skyModel = SkyModel.get_sky_model_from_fits(
                    # #     fits_file=fits_path,
                    # #     prefix_mapping=prefix_mapping,
                    # #     unit_mapping=unit_mapping,
                    # #     units_sources=units_sources,
                    # #     min_freq=None,
                    # #     max_freq=None,
                    # #     encoded_freq=None,
                    # #     memmap=False,
                    # # )
                    skyModel = create_sky_model_from_file(args.fits)
                    sources = skyModel.sources
                    print(f"Loaded {len(sources)} sources from {args.fits}")
                    
                else:
                    if args.json is not None:
                        # printlog (log_file, f"Loading sources from JSON file {args.json}")
                        # if (not os.path.isabs(args.json)):
                        #     fjson = os.path.join(os.path.dirname(__file__), args.json)
                        # else:
                        #     fjson = args.json
                        # with open(fjson, 'r') as f:
                        #     sources_data = json.load(f)
                        # sources = []
                        # for source_data in sources_data:
                        #     source = Source.from_json(source_data)
                        #     if args.scale_I is not None:
                        #         source.I *= args.scale_I
                        #     if source.ref_freq == 0:
                        #         source.ref_freq = args.ref_freq[0] * u.Hz if args.ref_freq is not None else frequency.to(u.Hz)
                        #     sources.append(source)
                        # if len(sources) == 0:
                        #     raise ValueError("No sources found in JSON file")
                        # source_ref = sources[0]
                        skyModel = create_sky_model_from_file(args.json)
                        sources = skyModel.sources
    

                    else:
                        source_ref = Source.from_name('HCG16')
                        N_sources = len(args.I)
                        intensities = args.I * u.Jy
                        sources_names = ['Random_{i+1}' for i in range(N_sources)]
                        sources = []
                        intensities = args.I
                        if intensities is not None:
                            if len(intensities) != len(sources_names):
                                for i in range(len(sources_names) - len(intensities)):
                                    intensities.append(0)
                            intensities = [intensity * u.Jy for intensity in intensities]
                        else:
                            intensities = [np.random.randint(1,10) * u.Jy for _ in range(len(sources_names))]

                        for idx,source_name in enumerate(sources_names):
                            try:
                                if idx == 0:
                                    source= source_ref
                                    source.I = intensities[idx]
                                else:
                                    x_coord = np.random.uniform(-fov.value/2, fov.value/2) * 0.8 * u.rad
                                    y_coord = np.random.uniform(-fov.value/2, fov.value/2) * 0.8 * u.rad
                                    source = Source(source.ra + x_coord, source.dec + y_coord, intensities[idx],)
                                
                                sources.append(source)
                            except Exception as e:
                                printlog(log_file, show_exc(e))
                    skyModel = SkyModel()
                    sources_list = []
                    for source in sources:
                        sources_list.append(source.to_sky_model(reduced_form=True))

                    sources_list = np.array(sources_list)            
                    skyModel.add_point_sources(sources_list)

                if args.json_fg is not None:
                    sources_fg = []
                    printlog (log_file, f"Loading foreground sources from JSON file {args.json_fg}")
                    if (not os.path.isabs(args.json_fg)):
                        fjson = os.path.join(os.path.dirname(__file__), args.json_fg)
                    else:
                        fjson = args.json_fg
                    with open(fjson, 'r') as f:
                        sources_data = json.load(f)
                    for source_data in sources_data:
                        source = Source.from_json(source_data)
                        sources_fg.append(source)
                    if len(sources_fg) == 0:
                        raise ValueError("No foreground sources found in JSON file")
                    
                    skyModel_fg = SkyModel()
                    sources_list_fg = []
                    for source in sources_fg:
                        sources_list_fg.append(source.to_sky_model(reduced_form=True))
                    sources_list_fg = np.array(sources_list_fg)
                    skyModel_fg.add_point_sources(sources_list_fg)

            if args.center is not None:
                try:
                    coords_str = args.center.replace(","," ").replace(":"," ")
                    center = SkyCoord(coords_str, unit=(u.hourangle, u.deg))    
                    # skyModel.set_center(center)
                except Exception as e:
                    print(show_exc(e))
                    center = skyModel.get_center()
            else:
                center = skyModel.get_center()


            source_ref = Source(center.ra, center.dec, 1, 0)

            if (args.delta_freq == 0) and (args.n_channels == 0):
                printlog (log_file, "Error: You must provide at least one of the two parameters: number of channels, delta frequency")
                sys.exit(1)
            if (args.delta_freq != 0) and (args.n_channels != 0):
                printlog (log_file, "Error: You must provide only one of the two parameters: number of channels, delta frequency")
                sys.exit(1)

            if (args.bandwidth == 0):
                n_channels = args.n_channels
                bandwidth = args.delta_freq * n_channels * u.MHz
                delta_freq = args.delta_freq * u.MHz
            else:
                bandwidth = args.bandwidth * u.MHz
                if (args.n_channels == 0):
                    n_channels = int(args.bandwidth / args.delta_freq)
                    delta_freq = args.delta_freq * u.MHz
                else:
                    n_channels = args.n_channels
                    delta_freq = bandwidth / n_channels

            start_freq = args.freq * u.MHz - n_channels * delta_freq/ 2 

                


            printlog (log_file, "Starting simulation with params:")
            printlog( log_file, f"\tCenter: {center.to_string('hmsdms')}")
            printlog (log_file, f"\tSources: {len(sources)}")
            printlog (log_file, f"\tPrefix: {prefix}")
            printlog (log_file, f"\tTelescope: {args.telescope}")
            printlog (log_file, f"\t\tDiameter: {get_diameter(args.telescope.upper()).value} meters")
            printlog (log_file, f"\t\tFOV: {fov.to(u.deg).value} degrees")
            printlog (log_file, f"\tFrequency: {args.freq} MHz")
            printlog (log_file, f"\tBandwidth: {bandwidth.to(u.MHz).value} MHz")
            printlog (log_file, f"\tNumber of channels: {n_channels}")
            printlog (log_file, f"\tDelta frequency: {delta_freq.value} MHz")
            printlog (log_file, f"\tObservation time: {args.seconds} seconds")
            printlog (log_file, f"\tPixels: {args.pixels}")
            printlog (log_file, f"\tIntensity: {args.I}")
            printlog (log_file, f"\tQ: {args.Q}")
            printlog (log_file, f"\tU: {args.U}")
            printlog (log_file, f"\tV: {args.V}")
            printlog (log_file, f"\tRef freq: {args.ref_freq}")
            printlog (log_file, f"\tCleaning: {args.cleaning}")
            if args.cleaning:
                printlog (log_file, f"\t\tNiter: {args.niter}")
            printlog (log_file, f"\tRMS: {args.rms}")
            if args.rms:
                printlog (log_file, f"\t\tRMS value: {args.rms_value} Jy")
                printlog (log_file, f"\t\tRMS sigma: {args.rms_sigma}")                

            # Launch the scripts monitor.py (not blocking) with pid of the current process
            # monitor_script = os.path.join(os.path.dirname(__file__), 'test_monitor.py')
            # if os.path.exists(monitor_script):
            #     printlog (log_file, f"Launching monitor script {monitor_script} with PID {os.getpid()}")
            #     monitor_proc = subprocess.Popen([sys.executable, monitor_script, str(os.getpid()),f"--csv={prefix}_monitor.log"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, close_fds=True)



            observation_time = source_ref.get_best_observation_time(telescope=telescope)
            number_of_timesteps = max(1,int(args.seconds / 7.997))
            number_of_channels = n_channels

            observation = Observation(
                start_frequency_hz=start_freq.to(u.Hz).value,
                start_date_and_time=observation_time,
                frequency_increment_hz=delta_freq.to(u.Hz).value,
                length = timedelta(seconds=args.seconds),
                number_of_time_steps=number_of_timesteps,
                number_of_channels=n_channels,
                phase_centre_ra_deg = center.ra.to(u.deg).value,
                phase_centre_dec_deg = center.dec.to(u.deg).value,
            )

            root_path = work_dir
            visibility_path = os.path.join(root_path, f'{prefix}_visibilities.MS')
            if os.path.exists(visibility_path):
                if args.overwrite:
                    printlog (log_file, f"Visibility file {visibility_path} already exists. Overwriting it")
                    shutil.rmtree(visibility_path)
                else:
                # Ask por confirmation to overwrite
                    print (log_file, f"Visibility file {visibility_path} already exists. Do you want to overwrite it? (y/n)")
                    answer = input()
                    if answer.lower() != 'y':
                        printlog (log_file, f"No overwrite, then exiting")
                        sys.exit(0)
                    else:
                        printlog (log_file, f"Overwriting visibility file {visibility_path}")
                        shutil.rmtree(visibility_path)
            printlog (log_file, f"Creating visibility file {visibility_path}")

            
            simulation_params = {
                'channel_bandwidth_hz': delta_freq.to(u.Hz).value,
                #'station_type': "Aperture array",
                'station_type': "Gaussian beam",
                'gauss_beam_fwhm_deg': fov.to(u.deg).value,
                'gauss_ref_freq_hz': frequency.to(u.Hz).value,
                'use_gpus': False,
            } 

            if args.rms:
                simulation_params['noise_enable'] = True
                
                # simulation_params['noise_rms_end'] = args.rms_value
                # simulation_params['noise_freq'] = "Observation settings"
                # simulation_params['noise_rms_start'] = args.rms_start
                # simulation_params['noise_rms_end'] = args.rms_end
                simulation_params['noise_freq'] = "Telescope model"
                simulation_params['noise_rms'] = "Telescope model"

            simulation = InterferometerSimulation(**simulation_params)
            # skyModel.set_center(center)
            simulation.run_simulation(telescope=telescope, observation=observation, sky=skyModel, visibility_path=visibility_path, backend=backend)

            printlog (log_file, f"Visibilities saved in {visibility_path}")
            printlog (log_file, "Recovering visibilities")
            visibilities = Visibility(visibility_path)
            imaging_cellsize = fov / int(args.pixels)
            if not args.cleaning :
                config = OskarDirtyImagerConfig(
                    imaging_npixel=args.pixels,
                    imaging_cellsize= imaging_cellsize.to(u.rad).value,
                    combine_across_frequencies=True,
                    imaging_phase_centre=source_ref.coords())
                imager = OskarDirtyImager(config=config)
                dirty_image = imager.create_dirty_image(visibilities)
                dirty_png_path = os.path.join(root_path, f'{prefix}_dirty.png')
                dirty_image.plot(title=f"Dirty image {backend.name} ({telescope.name.upper()})", wcs_enabled=True, xlabel='RA', ylabel='DEC')

                dirty_image.plot(title=f"Dirty image {backend.name} ({telescope.name.upper()})", filename=dirty_png_path, wcs_enabled=True, xlabel='RA', ylabel='DEC')
                printlog (log_file, f"Dirty image (PNG) saved in {dirty_png_path}")

                dirty_fits_path = os.path.join(root_path, f'{prefix}_dirty.fits')
                dirty_image.write_to_file(dirty_fits_path, overwrite=True)
                printlog (log_file, f"Dirty image (FITS) saved in {dirty_fits_path}")

            if args.cleaning:
                printlog (log_file, "Cleaning not supported for OSKAR, using WSCLEAN")
                path_fits = os.path.join(root_path, f"{prefix}_cleaned.fits")
                print (f"PATH Fits: {path_fits}")
                threshold_settings = "--auto-threshold 0.3"
                # if (args.rms) and (args.rms_value > 0):
                #     threshold_settings = f" --threshold {args.rms_value * args.rms_sigma}"
                custom_command = f"wsclean -weight briggs {args.robust} -multiscale -size {args.pixels} {args.pixels} -scale {imaging_cellsize.to(u.deg).value}deg -niter {args.niter} -mgain 0.8 {threshold_settings} -auto-mask 3 -channels-out 8 -join-channels {visibility_path}"
                printlog (log_file, f"Running custom wsclean command: {custom_command}")
                cleaned = iaa_create_image_custom_command(custom_command, path_fits)
                # Remove the temporary files created by WSClean
                tmp_files = glob.glob(os.path.join(tmp_dir, "wsclean-00*.fits"))
                for tmp_file in tmp_files:
                    os.remove(tmp_file)
              
                gamma = 0.3
                wsclean_files = glob.glob('wsclean-*.fits')
                seconds = args.seconds
                for img_path in wsclean_files:
                    img = Image(path=img_path)
                    img.plot(title=f"Cleaned image (WSCLEAN) {backend.name.upper()} ({telescope.name.upper()})", filename=f"{prefix}_{img_path.replace('fits','png')}", wcs_enabled=True, xlabel='RA', ylabel='DEC', norm=PowerNorm(gamma))
                    new_path = f'{prefix}_bw{bandwidth.to(u.MHz).value:.0f}_ch{n_channels}_fr{frequency.to(u.MHz).value:.0f}_sec{seconds}{img_path}'
                    new_path = new_path.replace('wsclean-', '')
                    shutil.move(img_path, new_path)


            # printlog(log_file, "=== Foreground simulation ===")
            # visibility_path = os.path.join(root_path, f'{prefix}_fg_visibilities.MS')
            # simulation.simulate_foreground_vis(telescope=telescope, observation=observation, sky=skyModel, visibility_path=visibility_path, backend=backend)

            # printlog (log_file, f"Visibilities saved in {visibility_path}")
            # printlog (log_file, "Recovering visibilities")
            # visibilities = Visibility(visibility_path)
            # imaging_cellsize = fov / int(args.pixels)
            # if not args.cleaning :
            #     config = OskarDirtyImagerConfig(
            #         imaging_npixel=args.pixels,
            #         imaging_cellsize= imaging_cellsize.to(u.rad).value,
            #         combine_across_frequencies=True,
            #         imaging_phase_centre=source_ref.coords())
            #     imager = OskarDirtyImager(config=config)
            #     dirty_image = imager.create_dirty_image(visibilities)
            #     dirty_png_path = os.path.join(root_path, f'{prefix}_dirty_fg.png')
            #     dirty_image.plot(title=f"Dirty image {backend.name} ({telescope.name.upper()})", wcs_enabled=True, xlabel='RA', ylabel='DEC')

            #     dirty_image.plot(title=f"Dirty image {backend.name} ({telescope.name.upper()})", filename=dirty_png_path, wcs_enabled=True, xlabel='RA', ylabel='DEC')
            #     printlog (log_file, f"Dirty image (PNG) saved in {dirty_png_path}")

            #     dirty_fits_path = os.path.join(root_path, f'{prefix}_dirty_fg.fits')
            #     dirty_image.write_to_file(dirty_fits_path, overwrite=True)
            #     printlog (log_file, f"Dirty image (FITS) saved in {dirty_fits_path}")

            # if args.cleaning:
            #     printlog (log_file, "Cleaning not supported for OSKAR, using WSCLEAN")
            #     path_fits = os.path.join(root_path, f"{prefix}_cleaned_fg.fits")
            #     threshold_settings = "--auto-threshold 0.3"
            #     # if (args.rms) and (args.rms_value > 0):
            #     #     threshold_settings = f" --threshold {args.rms_value * args.rms_sigma}"
            #     custom_command = f"wsclean -weight briggs {args.robust} -multiscale -size {args.pixels} {args.pixels} -scale {imaging_cellsize.to(u.deg).value}deg -niter {args.niter} -mgain 0.8 {threshold_settings} -auto-mask 3 -channels-out 8 -join-channels {visibility_path}"
            #     printlog (log_file, f"Running custom wsclean command: {custom_command}")
            #     cleaned = iaa_create_image_custom_command(custom_command, path_fits)
            #     # Remove the temporary files created by WSClean
            #     tmp_files = glob.glob(os.path.join(tmp_dir, "wsclean-00*.fits"))
            #     for tmp_file in tmp_files:
            #         os.remove(tmp_file)
              
            #     gamma = 0.3
            #     wsclean_files = glob.glob('wsclean-*.fits')
            #     seconds = args.seconds
            #     for img_path in wsclean_files:
            #         img = Image(path=img_path)
            #         img.plot(title=f"Cleaned image (WSCLEAN) {backend.name.upper()} ({telescope.name.upper()})", filename=f"{prefix}_{img_path.replace('fits','png')}", wcs_enabled=True, xlabel='RA', ylabel='DEC', norm=PowerNorm(gamma))
            #         new_path = f'{prefix}_bw{bandwidth.to(u.MHz).value:.0f}_ch{n_channels}_fr{frequency.to(u.MHz).value:.0f}_sec{seconds}{img_path}_fg'
            #         new_path = new_path.replace('wsclean-', '')
            #         shutil.move(img_path, new_path)         
            printlog (log_file, "Time ellapsed: ", time.time() - t0)




        except Exception as e:
            printlog (log_file, show_exc(e)) 
            printlog (log_file, "Time ellapsed: ", time.time() - t0)

            sys.exit(1)
    except Exception as e:
        print(show_exc(e))
        print ("Time ellapsed: ", time.time() - t0)
        sys.exit(1)
