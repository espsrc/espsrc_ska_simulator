"""cli.py — command-line entry point (adapted from synthsim.py)"""

# suppress xarray accessor collision between rascil and ska-sdp-datamodels
import warnings

warnings.filterwarnings("ignore", message="registration of accessor")

import argparse
import sys
from typing import List, Optional

from .config import ImgConfig, ObsConfig, SimConfig
from .pipeline import run


def _list_str_to_floats(s: Optional[str]) -> Optional[List[float]]:
    if s is None:
        return None
    return [float(x) for x in s.split(",")]


def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(
        description="SKA simulator — sky -> visibilities -> image"
    )
    p.add_argument(
        "--model", type=str, default=None, help="Sky model file (pickle/FITS/JSON)"
    )
    p.add_argument(
        "--fits", dest="fits_file", type=str, default=None, help="FITS catalogue file"
    )
    p.add_argument(
        "--json", dest="json_file", type=str, default=None, help="JSON catalogue file"
    )
    p.add_argument("--json-fg", type=str, default=None, help="JSON foreground sources")
    p.add_argument(
        "--center",
        type=str,
        default=None,
        help="Field centre, e.g. '10h01m35.1s 2d41m41s'",
    )
    p.add_argument("--prefix", type=str, default=None, help="Output prefix")
    p.add_argument("--telescope", type=str, default="SKA1MID", help="Telescope name")
    p.add_argument(
        "--telescope-version", type=str, default=None, help="Explicit telescope version"
    )
    p.add_argument(
        "--I", type=float, nargs="+", default=[10.0], help="Source intensities (Jy)"
    )
    p.add_argument("--Q", type=float, default=None, help="Stokes Q")
    p.add_argument("--U", type=float, default=None, help="Stokes U")
    p.add_argument("--V", type=float, default=None, help="Stokes V")
    p.add_argument(
        "--ref-freq",
        type=float,
        nargs="+",
        default=None,
        help="Reference freq(s) in Hz",
    )
    p.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing outputs"
    )
    p.add_argument("--freq", type=float, default=700.0, help="Centre frequency (MHz)")
    p.add_argument("--bandwidth", type=float, default=None, help="Bandwidth (MHz)")
    p.add_argument("--n-channels", type=int, default=None, help="Number of channels")
    p.add_argument("--delta-freq", type=float, default=None, help="Channel width (MHz)")
    p.add_argument("--seconds", type=int, default=10, help="Observation time (s)")
    p.add_argument(
        "--cleaning", action="store_true", help="Use WSClean instead of OSKAR dirty"
    )
    p.add_argument("--pixels", type=int, default=512, help="Image size in pixels")
    p.add_argument(
        "--imaging-niter", type=int, default=1000, help="OSKAR imaging iterations"
    )
    p.add_argument("--fov", type=str, default=None, help="Field of view, e.g. '1deg'")
    p.add_argument("--robust", type=float, default=0.0, help="Briggs robustness")
    p.add_argument("--rms", action="store_true", default=False, help="Enable noise")
    p.add_argument("--rms-value", type=float, default=0.0, help="RMS noise level (Jy)")
    p.add_argument("--rms-sigma", type=float, default=3.0, help="RMS sigma multiplier")
    p.add_argument("--niter", type=int, default=5000, help="WSClean iterations")
    p.add_argument("--scale-I", type=float, default=1.0, help="Scale Stokes I")
    p.add_argument(
        "--catalogue", type=int, default=0, help="1=MIGHTEE, 2=GLEAM, 3=SKAMid"
    )
    p.add_argument(
        "--column-mapping",
        type=str,
        default="0,1,2,3,4,5,6,7,8,9,10,11,12",
        help="FITS column indices: id,ra,dec,I,Q,U,V,spec_index,ref_freq,rm,major,minor,pa",
    )
    p.add_argument(
        "--show-telescopes",
        action="store_true",
        help="List available telescopes and exit",
    )

    args = p.parse_args(argv)

    if args.show_telescopes:
        from typing import get_args

        from karabo.simulation.telescope import (
            OSKARTelescopesWithoutVersionType,
            OSKARTelescopesWithVersionType,
        )

        choices = get_args(OSKARTelescopesWithVersionType) + get_args(
            OSKARTelescopesWithoutVersionType
        )
        for t in choices:
            print(t)
        sys.exit(0)

    # map CLI args to SimConfig
    obs = ObsConfig(
        freq_mhz=args.freq,
        bandwidth_mhz=args.bandwidth,
        n_channels=args.n_channels,
        delta_freq_mhz=args.delta_freq,
        seconds=args.seconds,
    )

    fov_deg = None
    if args.fov:
        try:
            import astropy.units as u

            fov_deg = u.Quantity(args.fov).to(u.deg).value
        except Exception:
            fov_deg = float(args.fov)

    img = ImgConfig(
        pixels=args.pixels,
        fov_deg=fov_deg,
        imaging_niter=args.imaging_niter,
        robust=args.robust,
        algorithm="wsclean_clean" if args.cleaning else "oskar_dirty",
    )

    sky_file = args.model or args.fits_file or args.json_file

    config = SimConfig(
        telescope=args.telescope,
        telescope_version=args.telescope_version,
        sky_file=sky_file,
        sky_format="auto",
        catalogue=args.catalogue,
        column_mapping=args.column_mapping,
        scale_I=args.scale_I,
        I=args.I,
        Q=args.Q,
        U=args.U,
        V=args.V,
        ref_freq_hz=args.ref_freq,
        json_fg=args.json_fg,
        center=args.center,
        rms=args.rms,
        rms_value=args.rms_value,
        rms_sigma=args.rms_sigma,
        niter=args.niter,
        observation=obs,
        imaging=img,
        output_prefix=args.prefix,
        overwrite=args.overwrite,
        cleaning=args.cleaning,
    )

    run(config)


if __name__ == "__main__":
    main()
