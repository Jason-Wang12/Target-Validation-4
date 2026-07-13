"""rcvtc pipeline driver.

Usage:
    python run_pipeline.py --config config/example_config.yaml
    python run_pipeline.py --config config/example_config.yaml --stages 1,2
    python run_pipeline.py --config config/example_config.yaml --resume

Stage list:
    1 = variant catalog + associations
    2 = functional characterization
    3 = functional-bin allelic series
    4 = published allelic-series retrieval
    5 = FinnGen replication
    6 = constraint + mechanism verdict
    7 = curated evidence + literature
    8 = expression + tissue
    9 = synthesis (cards + report)
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

from rcvtc.io_config import load_config, dump_run_config
from rcvtc.utils.logging import setup_logging, stage_logger, Checkpoint
from rcvtc.utils.ids import resolve_many

STAGE_MAP = {}   # populated below as stage modules land


def register(n: int):
    def deco(fn):
        STAGE_MAP[n] = fn
        return fn
    return deco


@register(1)
def stage1(cfg, ctx, ckpt):
    from rcvtc.stage1_variants import run_stage1
    return run_stage1(cfg, ctx, ckpt)


@register(2)
def stage2(cfg, ctx, ckpt):
    from rcvtc.stage2_functional import run_stage2
    return run_stage2(cfg, ctx, ckpt)


@register(3)
def stage3(cfg, ctx, ckpt):
    from rcvtc.stage3_allelic_series import run_stage3
    return run_stage3(cfg, ctx, ckpt)


@register(4)
def stage4(cfg, ctx, ckpt):
    from rcvtc.stage4_evidence import run_stage4
    return run_stage4(cfg, ctx, ckpt)


@register(5)
def stage5(cfg, ctx, ckpt):
    from rcvtc.stage5_finngen import run_stage5
    return run_stage5(cfg, ctx, ckpt)


@register(6)
def stage6(cfg, ctx, ckpt):
    from rcvtc.stage6_constraint import run_stage6
    return run_stage6(cfg, ctx, ckpt)


@register(7)
def stage7(cfg, ctx, ckpt):
    from rcvtc.stage7_curation import run_stage7
    return run_stage7(cfg, ctx, ckpt)


@register(8)
def stage8(cfg, ctx, ckpt):
    from rcvtc.stage8_expression import run_stage8
    return run_stage8(cfg, ctx, ckpt)


@register(9)
def stage9(cfg, ctx, ckpt):
    from rcvtc.stage9_synthesis import run_stage9
    return run_stage9(cfg, ctx, ckpt)


def _parse_stages(spec: str | None) -> list[int]:
    if not spec:
        return sorted(STAGE_MAP.keys())
    out: set[int] = set()
    for token in spec.split(","):
        token = token.strip()
        if "-" in token:
            a, b = token.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(token))
    return sorted(s for s in out if s in STAGE_MAP)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="Path to YAML config")
    p.add_argument("--stages", default=None, help="Comma-separated stage numbers or ranges, e.g. 1,2,3 or 1-4")
    p.add_argument("--resume", action="store_true", help="Skip stages with existing .done marker")
    p.add_argument("--force", action="store_true", help="Rerun stages even if .done exists")
    args = p.parse_args()

    cfg = load_config(args.config)
    outdir = Path(cfg["outputs_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    log = setup_logging(outdir)
    log.info("rcvtc pipeline start | run_id=%s | outputs=%s | config_hash=%s",
             cfg["run_id"], outdir, cfg["_config_hash"])

    dump_run_config(cfg, outdir / "run_config.snapshot.yaml")

    # Resolve gene IDs up-front so every stage can share them
    ctx: dict = {"gene_ids": resolve_many(cfg["genes"])}
    for sym, g in ctx["gene_ids"].items():
        if g is None:
            log.error("Failed to resolve gene: %s (halting)", sym)
            return 2
        log.info("gene resolved | %s | %s | %s", g.symbol, g.ensembl_id, g.uniprot_id)

    stages = _parse_stages(args.stages)
    log.info("stages to run: %s", stages)

    ckpt = Checkpoint(outdir)
    t0 = time.time()
    for n in stages:
        stage_name = f"stage{n}"
        if ckpt.exists(stage_name) and args.resume and not args.force:
            log.info("stage %d: already done (resume) — skipping", n)
            continue
        if args.force:
            ckpt.clear(stage_name)
        t_stage = time.time()
        log.info("stage %d: start", n)
        try:
            meta = STAGE_MAP[n](cfg, ctx, ckpt) or {}
            ckpt.mark(stage_name, meta={"wall_s": round(time.time() - t_stage, 2), **meta})
            log.info("stage %d: done in %.1fs", n, time.time() - t_stage)
        except Exception:
            log.exception("stage %d: FAILED", n)
            return 3

    log.info("pipeline done | total_wall=%.1fs", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
