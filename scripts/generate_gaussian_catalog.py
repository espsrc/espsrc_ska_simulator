import argparse
import json
import os

import numpy as np


def generate_catalog(seed: int = 42, n_sources: int = 100):
    rng = np.random.default_rng(seed)
    center_ra = 150.0
    center_dec = 2.0

    sources = []
    for i in range(n_sources):
        # 1 sq degree spread
        ra = center_ra + rng.uniform(-0.5, 0.5)
        dec = center_dec + rng.uniform(-0.5, 0.5)

        major_axis = sample_major_axis_arcsec(rng)
        axis_ratio = rng.uniform(0.25, 1.0)
        minor_axis = max(0.1, major_axis * axis_ratio)
        pa = rng.uniform(0, 180)
        spec_index = float(rng.normal(loc=-0.5, scale=0.2))

        # Flux density
        flux_i = sample_flux_jy(rng, i)

        src = {
            "ra": ra,
            "dec": dec,
            "I": flux_i,
            "Q": 0.0,
            "U": 0.0,
            "V": 0.0,
            "ref_freq": 1400000000.0,
            "spec_index": spec_index,
            "rot_meas": 0.0,
            "major_axis": major_axis,
            "minor_axis": minor_axis,
            "pa": pa,
            "true_redshift": 0.0,
            "obs_redshift": 0.0,
        }
        sources.append(src)

    os.makedirs("demo_output", exist_ok=True)
    out_path = "demo_output/reference_gaussian_catalog.json"
    with open(out_path, "w") as f:
        json.dump(sources, f, indent=2)
    ds9_path = "demo_output/reference_gaussian_catalog.reg"
    write_ds9_regions(sources, ds9_path)
    print(f"Generated {n_sources} Gaussian sources with seed={seed} in {out_path}")
    print(f"Generated DS9 region file in {ds9_path}")


def sample_major_axis_arcsec(rng):
    """Draw a broad source-size distribution in arcsec."""
    branch = rng.random()
    if branch < 0.12:
        return float(rng.uniform(0.1, 1.0))
    if branch < 0.62:
        return float(np.clip(rng.lognormal(mean=np.log(3.0), sigma=0.45), 1.0, 15.0))
    if branch < 0.85:
        return float(rng.uniform(15.0, 60.0))
    return float(rng.uniform(60.0, 300.0))


def sample_flux_jy(rng, index):
    """Draw a demo flux distribution centered near 1 mJy."""
    if index < 2:
        return float(rng.uniform(0.06, 0.10))
    branch = rng.random()
    if branch < 0.18:
        return float(10 ** rng.uniform(np.log10(1.0e-5), np.log10(2.0e-5)))
    if branch < 0.92:
        return float(10 ** rng.normal(np.log10(1.0e-3), 0.35))
    return float(10 ** rng.uniform(np.log10(5.0e-3), np.log10(2.0e-2)))


def write_ds9_regions(sources, path):
    """Write source ellipses to a DS9 FK5 region file."""
    lines = [
        "# Region file format: DS9 version 4.1",
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" '
        "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1",
        "fk5",
    ]
    for idx, src in enumerate(sources, start=1):
        ds9_angle = ds9_angle_from_catalog_pa(src["pa"])
        lines.append(
            "ellipse("
            f"{src['ra']:.10f},{src['dec']:.10f},"
            f'{src["major_axis"]:.6f}",{src["minor_axis"]:.6f}",'
            f"{ds9_angle:.6f}"
            f") # text={{src_{idx:03d} I={src['I']:.4g} Jy alpha={src['spec_index']:.3f}}}"
        )
        lines.append(
            f"point({src['ra']:.10f},{src['dec']:.10f}) "
            f"# point=cross 8 color=cyan text={{src_{idx:03d}}}"
        )
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def ds9_angle_from_catalog_pa(pa_deg):
    """Convert catalog PA east of north to DS9 display angle from +x."""
    return (90.0 - float(pa_deg)) % 180.0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a deterministic reference sky model catalog."
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Seed for random number generator"
    )
    parser.add_argument(
        "--n-sources", type=int, default=100, help="Number of sources to generate"
    )
    args = parser.parse_args()
    generate_catalog(seed=args.seed, n_sources=args.n_sources)
