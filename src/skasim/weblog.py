"""weblog.py — render weblog.html from a RunManifest via Jinja2."""

from __future__ import annotations

import importlib.resources
import json
from pathlib import Path
import base64
from jinja2 import Environment, PackageLoader, select_autoescape

from .manifest import RunManifest


def _humanize_seconds(total_s: float) -> str:
    """convert seconds to '2m 34s' or '1h 2m 3s'."""
    if total_s < 60:
        return f"{total_s:.1f}s"
    minutes, seconds = divmod(int(total_s), 60)
    if minutes < 60:
        return f"{minutes}m {seconds}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m {seconds}s"


def _find_image_outputs(manifest: RunManifest, work_dir: Path) -> list[dict]:
    """extract image paths, captions, and base64 data for inline embedding."""
    results: list[dict] = []
    for o in manifest.outputs:
        if not o.endswith((".png", ".jpg", ".jpeg")):
            continue
        fpath = work_dir / o
        caption = _infer_image_caption(Path(o).name)
        data_b64 = _file_to_base64_data_uri(fpath)
        results.append({"path": o, "caption": caption, "data": data_b64})
    return results


def _file_to_base64_data_uri(fpath: Path) -> str:
    """read a binary image and return a base64 data URI string."""
    data = fpath.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    # infer mime type from extension
    ext = fpath.suffix.lower()
    mime = "image/png" if ext == ".png" else "image/jpeg"
    return f"data:{mime};base64,{b64}"

def _infer_image_caption(filename: str) -> str:
    """return a human-readable caption from an output filename."""
    fname = Path(filename).stem
    lower = fname.lower()
    if "dirty" in lower:
        if "oskar" in lower or "-dirty" in lower:
            return "Dirty image (OSKAR)"
        return "Dirty image"
    if "psf" in lower:
        return "Point spread function"
    if "residual" in lower:
        return "Residual"
    if "model" in lower:
        return "Component model"
    if "clean" in lower or "image" in lower:
        return "Cleaned image (WSClean)"
    if "telescope" in lower:
        return "Telescope layout"
    return fname  # fallback to filename

def render_weblog(manifest: RunManifest, work_dir: Path) -> str:
    """generate weblog.html content and write it to disk. returns the html string."""
    env = Environment(
        loader=PackageLoader("skasim", "templates"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("weblog.html.j2")

    # compute summary stats
    total_elapsed = None
    if manifest.completed_at and manifest.started_at:
        total_elapsed = (manifest.completed_at - manifest.started_at).total_seconds()

    milestone_lookup = {m.name: m for m in manifest.milestones}
    phase_a_duration = None
    if "phase_a_started" in milestone_lookup and "phase_a_completed" in milestone_lookup:
        phase_a_duration = (
            milestone_lookup["phase_a_completed"].timestamp_utc
            - milestone_lookup["phase_a_started"].timestamp_utc
        ).total_seconds()

    phase_b_duration = None
    if "phase_b_started" in milestone_lookup and "phase_b_completed" in milestone_lookup:
        phase_b_duration = (
            milestone_lookup["phase_b_completed"].timestamp_utc
            - milestone_lookup["phase_b_started"].timestamp_utc
        ).total_seconds()

    html = template.render(
        manifest=manifest,
        total_elapsed=_humanize_seconds(total_elapsed) if total_elapsed else None,
        phase_a_duration=_humanize_seconds(phase_a_duration) if phase_a_duration else None,
        phase_b_duration=_humanize_seconds(phase_b_duration) if phase_b_duration else None,
        images=_find_image_outputs(manifest, work_dir),
        json_dump=lambda obj: json.dumps(obj, indent=2, default=str),
    )

    weblog_path = work_dir / "weblog.html"
    weblog_path.write_text(html, encoding="utf-8")
    return html
