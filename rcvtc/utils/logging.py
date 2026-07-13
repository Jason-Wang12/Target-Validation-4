"""Structured logging + checkpoint tracking."""
from __future__ import annotations
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Optional


def setup_logging(run_dir: str | Path, level: int = logging.INFO) -> logging.Logger:
    """Configure the root logger with console + file handlers.

    Log files are written to local /workspace (S3 FUSE doesn't permit append-mode
    open); after the run, driver code should copy the log into `run_dir/logs/`.
    """
    run_dir = Path(run_dir)
    log_dir_local = Path("/workspace/rcvtc_logs") / run_dir.name
    log_dir_local.mkdir(parents=True, exist_ok=True)

    fmt = "%(asctime)s | %(levelname)-7s | %(name)-25s | %(message)s"
    formatter = logging.Formatter(fmt, datefmt="%H:%M:%S")

    root = logging.getLogger("rcvtc")
    root.handlers.clear()
    root.setLevel(level)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    log_path = log_dir_local / "pipeline.log"
    fh = logging.FileHandler(log_path)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    root.info("log file: %s", log_path)
    return root


def stage_logger(stage_name: str) -> logging.Logger:
    return logging.getLogger(f"rcvtc.{stage_name}")


class Checkpoint:
    """Simple stage-level checkpoint tracking."""

    def __init__(self, run_dir: str | Path):
        self.run_dir = Path(run_dir)
        self.cp_dir = self.run_dir / "checkpoints"
        self.cp_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, stage: str) -> Path:
        return self.cp_dir / f"{stage}.done"

    def exists(self, stage: str) -> bool:
        return self._path(stage).exists()

    def mark(self, stage: str, meta: Optional[dict] = None) -> None:
        payload = {"stage": stage, "finished_at": time.time(), "meta": meta or {}}
        self._path(stage).write_text(json.dumps(payload, indent=2))

    def clear(self, stage: str) -> None:
        p = self._path(stage)
        if p.exists():
            p.unlink()
