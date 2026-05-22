"""manifest.py — Pydantic models for run manifest, milestone tracking, and pipeline context."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, Field

from .config import SimConfig
from .utils import init_logger


class Milestone(BaseModel):
    """single checkpoint in a simulation run."""

    name: str
    status: Literal["started", "completed", "failed"]
    timestamp_utc: datetime
    elapsed_s: Optional[float] = None
    details: dict = Field(default_factory=dict)


class RunManifest(BaseModel):
    """canonical machine-readable record of one simulation run."""

    run_id: str
    status: Literal["running", "completed", "failed"] = "running"
    started_at: datetime
    completed_at: Optional[datetime] = None
    config: SimConfig
    milestones: list[Milestone] = Field(default_factory=list)
    outputs: list[str] = Field(default_factory=list)
    errors: list[str] = Field(default_factory=list)

    def add_milestone(
        self,
        name: str,
        status: Literal["started", "completed", "failed"],
        elapsed_s: Optional[float] = None,
        details: Optional[dict] = None,
    ) -> Milestone:
        """append a milestone and return it."""
        ms = Milestone(
            name=name,
            status=status,
            timestamp_utc=datetime.utcnow(),
            elapsed_s=elapsed_s,
            details=details or {},
        )
        self.milestones.append(ms)
        return ms

    def mark_completed(self) -> None:
        """mark the run as completed."""
        self.status = "completed"
        self.completed_at = datetime.utcnow()

    def mark_failed(self, error: str) -> None:
        """mark the run as failed and record the error."""
        self.status = "failed"
        self.errors.append(error)

    def model_dump_json(self, **kwargs) -> str:
        """serialize to pretty-printed JSON string."""
        return json.dumps(self.model_dump(mode="json"), indent=2, default=str)


class RunContext(BaseModel):
    """passed through all pipeline functions; bundles config, paths, and manifest."""

    config: SimConfig
    work_dir: Path
    manifest: RunManifest
    visibility_path: Path
    log_path: Path
    manifest_path: Path
    weblog_path: Path
    sky_file_resolved: Optional[Path] = None

    def save_manifest(self) -> None:
        """write the current manifest state to disk (overwrites)."""
        self.manifest_path.write_text(self.manifest.model_dump_json(), encoding="utf-8")

    def add_milestone(self, *args, **kwargs) -> Milestone:
        """convenience: add milestone to manifest and persist to disk."""
        ms = self.manifest.add_milestone(*args, **kwargs)
        self.save_manifest()
        return ms


def create_run_context(config: SimConfig) -> RunContext:
    """create work_dir, init logger, build RunContext with empty manifest."""
    prefix = config.output_prefix or datetime.now().strftime("%Y%m%d_%H%M")
    prefix = f"{prefix}_{config.telescope.replace('-', '_')}"
    work_dir = Path(prefix).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    log_file = str(work_dir / f"{prefix}.log")
    init_logger(log_file)

    manifest = RunManifest(
        run_id=prefix,
        started_at=datetime.utcnow(),
        config=config,
    )

    ctx = RunContext(
        config=config,
        work_dir=work_dir,
        manifest=manifest,
        visibility_path=work_dir / "visibilities.MS",
        log_path=Path(log_file),
        manifest_path=work_dir / "run_manifest.json",
        weblog_path=work_dir / "weblog.html",
        sky_file_resolved=None,
    )

    if config.sky_file is not None:
        fpath = config.sky_file
        if not os.path.isabs(fpath):
            fpath = os.path.join(os.getcwd(), fpath)
        ctx.sky_file_resolved = Path(fpath).resolve()

    ctx.save_manifest()
    return ctx
