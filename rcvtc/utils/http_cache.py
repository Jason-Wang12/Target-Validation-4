"""Shared HTTP session with on-disk caching for public REST/GraphQL endpoints."""
from __future__ import annotations
import os
from pathlib import Path
from typing import Optional
import requests
import requests_cache


CACHE_DIR = Path(os.environ.get("RCVTC_CACHE", "/workspace/rcvtc_cache"))
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Default TTL 7 days for API responses; can be overridden per session.
_DEFAULT_TTL_SEC = 7 * 24 * 3600

_session: Optional[requests_cache.CachedSession] = None


def get_session(ttl_sec: int = _DEFAULT_TTL_SEC) -> requests_cache.CachedSession:
    global _session
    if _session is None:
        _session = requests_cache.CachedSession(
            cache_name=str(CACHE_DIR / "http_cache"),
            backend="sqlite",
            expire_after=ttl_sec,
            allowable_methods=("GET", "POST"),
            allowable_codes=(200,),
            stale_if_error=True,
        )
        # Standard headers
        _session.headers.update({
            "Accept": "application/json",
            "User-Agent": "rcvtc/0.1 (target-characterization; +https://github.com/phylo/rcvtc)",
        })
    return _session


def clear_cache() -> None:
    global _session
    if _session is not None:
        _session.cache.clear()
        _session = None
    p = CACHE_DIR / "http_cache.sqlite"
    if p.exists():
        p.unlink()
