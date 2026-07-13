"""Load, validate, and normalize the YAML config."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
from typing import Any, Optional
import yaml


class ConfigError(ValueError):
    """Raised when the config is malformed or references unknown resources."""


REQUIRED_TOP = {"run_id", "outputs_dir", "genes", "stage1", "stage2",
                "stage3", "stage4", "stage5", "stage6", "stage7", "stage8", "output"}


def load_config(path: str | Path) -> dict:
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")
    with p.open("r") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ConfigError(f"Config must be a top-level YAML mapping, got {type(cfg).__name__}")
    validate(cfg)
    cfg["_config_hash"] = _hash(cfg)
    cfg["_config_path"] = str(p.resolve())
    # Resolve outputs_dir with run_id substitution if present
    od = cfg["outputs_dir"]
    if "${run_id}" in od:
        cfg["outputs_dir"] = od.replace("${run_id}", cfg["run_id"])
    return cfg


def validate(cfg: dict) -> None:
    missing = REQUIRED_TOP - set(cfg.keys())
    if missing:
        raise ConfigError(f"Config missing required top-level keys: {sorted(missing)}")

    if not isinstance(cfg["genes"], list) or not cfg["genes"]:
        raise ConfigError("`genes` must be a non-empty list")
    if len(cfg["genes"]) > 20:
        raise ConfigError(f"`genes` list is capped at 20; got {len(cfg['genes'])}")

    for g in cfg["genes"]:
        if not isinstance(g, str) or not g.strip():
            raise ConfigError(f"Each gene must be a non-empty string; got {g!r}")

    # Stage 1 sub-checks
    s1 = cfg.get("stage1", {})
    af = s1.get("af_filter", {})
    if "max_popmax_af" not in af:
        raise ConfigError("stage1.af_filter.max_popmax_af is required")

    # Certainty enforcement is on by default
    out = cfg.get("output", {})
    if not out.get("certainty_labels_required_on_every_claim", True):
        raise ConfigError("output.certainty_labels_required_on_every_claim cannot be disabled")


def _hash(cfg: dict) -> str:
    """Stable hash of config (excluding runtime fields)."""
    cp = {k: v for k, v in cfg.items() if not k.startswith("_")}
    payload = json.dumps(cp, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def dump_run_config(cfg: dict, path: str | Path) -> None:
    """Snapshot the effective config to `path` for reproducibility.

    All runtime fields (leading `_`) are preserved. YAML round-trip.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # yaml.safe_dump won't handle arbitrary objects; make sure everything is
    # JSON-round-trip-safe first (converts dataclasses, pathlib, etc).
    safe = json.loads(json.dumps(cfg, default=str))
    with p.open("w") as f:
        yaml.safe_dump(safe, f, sort_keys=False)
