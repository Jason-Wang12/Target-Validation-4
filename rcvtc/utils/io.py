"""S3 FUSE-safe file I/O helpers.

S3-backed mounts (`/mnt/results/`, `/mnt/shared-workspace/`) do not permit some
POSIX operations that PyArrow / stdlib logging depend on. Workaround pattern:
    stage to /workspace/, then shutil.copy() to the S3 destination.

Formats we currently need this for:
    - Parquet (pyarrow)
    - Any other random-access binary write later (h5, sqlite, ...)
"""
from __future__ import annotations
import shutil
import tempfile
from pathlib import Path

import pandas as pd


def safe_write_parquet(df: pd.DataFrame, dest: str | Path) -> Path:
    """Write `df` to `dest` in a way that survives S3 FUSE.

    Writes to a temp file under /workspace/, then shutil.copy()'s into place.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    tmp_root = Path("/workspace/rcvtc_tmp")
    tmp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".parquet", dir=tmp_root, delete=False) as f:
        tmp_path = Path(f.name)

    try:
        df.to_parquet(tmp_path, index=False)
        shutil.copy(str(tmp_path), str(dest))
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
    return dest
