"""Create a flat spectral-index FITS map matching a Stokes-I model FITS file."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create an alpha FITS image with the same shape/WCS as a model FITS."
    )
    parser.add_argument("stokes_i", help="Input Stokes-I FITS image")
    parser.add_argument("alpha", help="Output spectral-index FITS image")
    parser.add_argument(
        "--value",
        type=float,
        default=-0.7,
        help="Flat spectral-index value to write",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output file",
    )
    args = parser.parse_args()

    stokes_i = Path(args.stokes_i)
    alpha = Path(args.alpha)
    alpha.parent.mkdir(parents=True, exist_ok=True)

    with fits.open(stokes_i) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()
    header["BUNIT"] = "1"
    header["COMMENT"] = "Flat spectral-index map generated for skasim testing."
    fits.writeto(
        alpha,
        np.full(data.shape, args.value, dtype=float),
        header,
        overwrite=args.overwrite,
    )
    print(f"Created {alpha} — flat alpha={args.value}, shape={data.shape}")


if __name__ == "__main__":
    main()
