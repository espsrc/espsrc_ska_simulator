#!/usr/bin/env python3
"""Quick shadems uv-coverage verification for demo visibility sets."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import time
from pathlib import Path


DEFAULT_PATTERNS = (
    "demo_outputs/tests*/visibilities.MS",
    "demo_outputs/test*/visibilities.MS",
    "demo_output/tests*/visibilities.MS",
    "demo_output/test*/visibilities.MS",
)


def slugify(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value.strip("_") or "visibility"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run shadems uv-coverage plots for demo visibilities and time each run."
    )
    parser.add_argument(
        "patterns",
        nargs="*",
        default=DEFAULT_PATTERNS,
        help="Glob(s) for Measurement Sets. Defaults cover demo_outputs/test* and demo_output/test*.",
    )
    parser.add_argument(
        "--output-dir",
        default="demo_output/shadems_uvcoverage_check",
        help="Directory for PNGs, logs, and timings CSV.",
    )
    parser.add_argument(
        "--shadems",
        default="shadems",
        help="shadems executable or path. Run this script inside the emcp conda env.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N matched Measurement Sets.",
    )
    parser.add_argument(
        "--row-chunk-size",
        default=None,
        help="Optional value passed to shadems --row-chunk-size.",
    )
    parser.add_argument(
        "--num-parallel",
        default=None,
        help="Optional value passed to shadems --num-parallel.",
    )
    parser.add_argument(
        "--canvas-size",
        type=int,
        default=600,
        help="Square shadems canvas size in pixels for U/V plots.",
    )
    return parser.parse_args()


def discover_measurement_sets(patterns: list[str]) -> list[Path]:
    matches: set[Path] = set()
    for pattern in patterns:
        for path in Path().glob(pattern):
            if path.is_dir():
                matches.add(path)
    return sorted(matches)


def test_name(ms: Path) -> str:
    return ms.parent.name


def directory_size(path: Path) -> int:
    total = 0
    for item in path.rglob("*"):
        if item.is_file():
            total += item.stat().st_size
    return total


def main() -> int:
    args = parse_args()

    if shutil.which(args.shadems) is None and not Path(args.shadems).exists():
        print(
            f"Could not find {args.shadems!r}. Activate the emcp conda env first, e.g. "
            "conda run -n emcp python scripts/verify_shadems_uvcoverage.py",
        )
        return 2

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Avoid slow or failing imports when user-level cache directories are not writable.
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(output_dir / ".mplconfig"))
    env.setdefault("NUMBA_CACHE_DIR", str(output_dir / ".numba_cache"))
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

    measurement_sets = discover_measurement_sets(list(args.patterns))
    if args.limit is not None:
        measurement_sets = measurement_sets[: args.limit]

    if not measurement_sets:
        print("No Measurement Sets matched these patterns:")
        for pattern in args.patterns:
            print(f"  {pattern}")
        return 1

    csv_path = output_dir / "uvcoverage_timings.csv"
    rows: list[dict[str, str]] = []

    for index, ms in enumerate(measurement_sets, start=1):
        name = test_name(ms)
        slug = slugify(name)
        png_name = f"{slug}-uvcoverage.png"
        log_path = output_dir / f"{slug}.log"
        ms_size_bytes = directory_size(ms)

        cmd = [
            args.shadems,
            str(ms),
            "-x",
            "u",
            "-y",
            "v",
            "--dir",
            str(output_dir),
            "--png",
            png_name,
            "--title",
            f"{name} uv coverage",
            "--xlabel",
            "u",
            "--ylabel",
            "v",
            "--xcanvas",
            str(args.canvas_size),
            "--ycanvas",
            str(args.canvas_size),
            "--spread-pix",
            "2",
            "--no-lim-save",
        ]
        if args.row_chunk_size:
            cmd.extend(["--row-chunk-size", args.row_chunk_size])
        if args.num_parallel:
            cmd.extend(["--num-parallel", args.num_parallel])

        print(f"[{index}/{len(measurement_sets)}] {name}: {ms}")
        started = time.perf_counter()
        completed = subprocess.run(
            cmd,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        status = "ok" if completed.returncode == 0 else "failed"
        returncode = str(completed.returncode)
        elapsed = time.perf_counter() - started
        elapsed_seconds = round(elapsed)
        log_path.write_text(completed.stdout, encoding="utf-8")

        png_path = output_dir / png_name
        row = {
            "test": name,
            "measurement_set": str(ms),
            "ms_size_bytes": str(ms_size_bytes),
            "ms_size_gib": f"{ms_size_bytes / 1024**3:.2f}",
            "seconds": str(elapsed_seconds),
            "status": status,
            "returncode": returncode,
            "png": str(png_path if png_path.exists() else ""),
            "log": str(log_path),
        }
        rows.append(row)
        print(f"  {row['status']} in {row['seconds']} s; MS size {row['ms_size_gib']} GiB")

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "test",
                "measurement_set",
                "ms_size_bytes",
                "ms_size_gib",
                "seconds",
                "status",
                "returncode",
                "png",
                "log",
            ),
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote timing summary to {csv_path}")
    return 0 if all(row["status"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
