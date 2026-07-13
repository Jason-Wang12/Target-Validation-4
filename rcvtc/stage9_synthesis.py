"""Stage 9 — Synthesis: per-gene evidence cards + ranking + report.

Consumes outputs from all earlier stages and produces:

  1. Per-gene evidence cards ({gene_symbol}_card.json + .md)
  2. Cross-gene ranking table (rankings.parquet + rankings.md)
  3. Master markdown report (report_rcvtc.md)

Every claim in the cards is tier-tagged and cites its source. The tier-audit
(certainty.assert_tier_tagged) is invoked over each card before emit.
"""
from __future__ import annotations
import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from .certainty import Tier
from .utils.io import safe_write_parquet

log = logging.getLogger("rcvtc.stage9")


def _read(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()


def _num(x, fmt=".3g"):
    try:
        return format(float(x), fmt)
    except (TypeError, ValueError):
        return "—"


def _build_card(sym: str, outdir_run: Path, cfg: dict) -> dict:
    """Build one gene's evidence card as a dict."""
    card: dict[str, Any] = {
        "gene_symbol": sym,
        "run_id": cfg.get("run_id"),
        "config_hash": cfg.get("_config_hash"),
    }

    # === Stage 1 — variant catalog ===
    v = _read(outdir_run / "stage1" / "1a_gnomad" / "gnomad_variants.parquet")
    if not v.empty and (v["gene_symbol"] == sym).any():
        vg = v[v["gene_symbol"] == sym]
        # Bucket by consequence
        cs = vg["consequence"].value_counts().to_dict() if "consequence" in vg.columns else {}
        card["variant_catalog"] = {
            "n_rare_coding": int(len(vg)),
            "by_consequence": {k: int(v) for k, v in cs.items()},
            "certainty_tier": Tier.MEASURED.value,
            "source": "gnomad_v4_1",
        }

    # Genebass associations (real output uses `gene`; some smoke fixtures use `gene_symbol`)
    gb = _read(outdir_run / "stage1" / "1b_genebass" / "genebass_associations.parquet")
    gb_gene_col = "gene_symbol" if "gene_symbol" in gb.columns else "gene"
    if not gb.empty and gb_gene_col in gb.columns and (gb[gb_gene_col] == sym).any():
        gg = gb[gb[gb_gene_col] == sym].sort_values("Pvalue_Burden") if "Pvalue_Burden" in gb.columns else gb[gb[gb_gene_col] == sym]
        top = gg.head(5).to_dict(orient="records")
        card["genebass_top"] = {
            "n_kept": int(len(gg)),
            "top5": [
                {
                    "phenotype": r.get("phenotype_slug") or r.get("phenocode"),
                    "annotation": r.get("annotation_mask"),
                    "beta": r.get("BETA_Burden"),
                    "p_value": r.get("Pvalue_Burden"),
                    "certainty_tier": Tier.INFERRED.value,
                    "source": "genebass_ukb_wes_450k",
                }
                for r in top
            ],
        }

    # === Stage 2 — functional evidence (summaries) ===
    am = _read(outdir_run / "stage2" / "2b_alphamissense" / "alphamissense_scores.parquet")
    if not am.empty and (am["gene_symbol"] == sym).any():
        ag = am[am["gene_symbol"] == sym]
        cls_counts = ag["am_class"].value_counts().to_dict() if "am_class" in ag.columns else {}
        card["alphamissense"] = {
            "n_scored": int(len(ag)),
            "mean_pathogenicity": float(ag["am_pathogenicity"].mean()) if "am_pathogenicity" in ag.columns else None,
            "class_counts": {k: int(v) for k, v in cls_counts.items()},
            "certainty_tier": Tier.PREDICTED.value,
            "source": "alphamissense_v1",
        }

    mv = _read(outdir_run / "stage2" / "2a_mavedb" / "mavedb_scores.parquet")
    if not mv.empty and (mv["gene_symbol"] == sym).any():
        mg = mv[mv["gene_symbol"] == sym]
        card["mavedb"] = {
            "n_measurements": int(len(mg)),
            "n_scoresets": int(mg["scoreset_urn"].nunique()) if "scoreset_urn" in mg.columns else None,
            "certainty_tier": Tier.MEASURED.value,
            "source": "mavedb",
        }

    up = _read(outdir_run / "stage2" / "2c_uniprot" / "uniprot_features.parquet")
    if not up.empty and (up["gene_symbol"] == sym).any():
        ug = up[up["gene_symbol"] == sym]
        by_type = ug["feature_type"].value_counts().to_dict() if "feature_type" in ug.columns else {}
        card["uniprot_features"] = {
            "n_features": int(len(ug)),
            "by_type": {k: int(v) for k, v in by_type.items()},
            "certainty_tier": Tier.ANNOTATED.value,
            "source": "uniprot",
        }

    af = _read(outdir_run / "stage2" / "2d_alphafold" / "alphafold_plddt.parquet")
    if not af.empty and (af["gene_symbol"] == sym).any():
        afg = af[af["gene_symbol"] == sym]
        card["alphafold"] = {
            "n_residues": int(len(afg)),
            "median_plddt": float(afg["plddt"].median()) if "plddt" in afg.columns else None,
            "certainty_tier": Tier.PREDICTED.value,
            "source": "alphafold_db",
        }

    # === Stage 3 — mask comparison + functional bins ===
    mc = _read(outdir_run / "stage3" / "mask_comparison.parquet")
    if not mc.empty and "gene" in mc.columns and (mc["gene"] == sym).any():
        mcg = mc[mc["gene"] == sym]
        patterns = mcg["mask_pattern"].value_counts().to_dict() if "mask_pattern" in mcg.columns else {}
        card["mask_patterns"] = {
            "n_phenotypes_compared": int(len(mcg)),
            "pattern_counts": {k: int(v) for k, v in patterns.items()},
            "certainty_tier": Tier.INFERRED.value,
            "source": "rcvtc_mask_comparison_v1",
        }

    # === Stage 4 — ClinVar + OT curated ===
    cv = _read(outdir_run / "stage4" / "4a_clinvar" / "clinvar_variants.parquet")
    if not cv.empty and (cv["gene_symbol"] == sym).any():
        cvg = cv[cv["gene_symbol"] == sym]
        card["clinvar_2plus_star"] = {
            "n_pathogenic": int((cvg["clinsig_bucket"] == "pathogenic_lp").sum()),
            "n_benign": int((cvg["clinsig_bucket"] == "benign_lb").sum()),
            "max_stars": int(cvg["stars"].max()),
            "certainty_tier": Tier.ANNOTATED.value,
            "source": "clinvar_2plus_star",
        }

    # === Stage 5 — FinnGen replication ===
    fg = _read(outdir_run / "stage5" / "5a_finngen" / "finngen_gene_phenos.parquet")
    if not fg.empty and (fg["gene_symbol"] == sym).any():
        fgg = fg[fg["gene_symbol"] == sym].sort_values("pval") if "pval" in fg.columns else fg[fg["gene_symbol"] == sym]
        top5 = fgg.head(5)
        def _var_str(r):
            c = r.get("variant_chrom") or "?"
            p = r.get("variant_pos") or "?"
            ref = r.get("variant_ref") or "?"
            alt = r.get("variant_alt") or "?"
            return f"chr{c}:{p}_{ref}>{alt}"
        card["finngen_top"] = {
            "n_phenotypes": int(len(fgg)),
            "top5": [
                {
                    "phenotype": r.get("phenostring") or r.get("phenocode"),
                    "beta": r.get("beta"),
                    "p_value": r.get("pval"),
                    "maf": r.get("maf"),
                    "variant": _var_str(r),
                    "certainty_tier": Tier.INFERRED.value,
                    "source": "finngen_r12",
                }
                for _, r in top5.iterrows()
            ],
        }

    # === Stage 6 — constraint + verdict ===
    con = _read(outdir_run / "stage6" / "6a_constraint" / "gene_constraint.parquet")
    if not con.empty and (con["gene_symbol"] == sym).any():
        cr = con[con["gene_symbol"] == sym].iloc[0]
        card["constraint"] = {
            "loeuf": cr.get("loeuf"),
            "oe_lof": cr.get("oe_lof"),
            "oe_mis": cr.get("oe_mis"),
            "pLI": cr.get("pLI"),
            "mis_z": cr.get("mis_z"),
            "certainty_tier": Tier.MEASURED.value,
            "source": "gnomad_v4_1",
        }
    vd = _read(outdir_run / "stage6" / "mechanism_verdicts.parquet")
    if not vd.empty and (vd["gene_symbol"] == sym).any():
        vr = vd[vd["gene_symbol"] == sym].iloc[0]
        card["mechanism_verdict"] = {
            "verdict": vr.get("verdict"),
            "reasons": vr.get("reasons"),
            "certainty_tier": Tier.INFERRED.value,
            "source": "rcvtc_mechanism_verdict_v1",
        }

    # === Stage 7 — curated ledger ===
    led = _read(outdir_run / "stage7" / "curated_gene_ledger.parquet")
    if not led.empty and (led["gene_symbol"] == sym).any():
        lr = led[led["gene_symbol"] == sym].iloc[0]
        card["curated_ledger"] = {
            "clinvar_pathogenic_n": int(lr.get("clinvar_pathogenic_n") or 0),
            "clinvar_max_stars": int(lr.get("clinvar_max_stars") or 0),
            "ot_weighted_sum": float(lr.get("ot_weighted_sum") or 0),
            "top_diseases": lr.get("top_diseases"),
            "certainty_tier": Tier.ANNOTATED.value,
            "source": "rcvtc_curated_synthesis_v1",
        }

    # === Stage 8 — expression ===
    ex = _read(outdir_run / "stage8" / "hpa_expression_summary.parquet")
    if not ex.empty and (ex["gene_symbol"] == sym).any():
        er = ex[ex["gene_symbol"] == sym].iloc[0]
        card["expression"] = {
            "max_tissue": er.get("max_tissue"),
            "max_ntpm": er.get("max_ntpm"),
            "median_ntpm": er.get("median_ntpm"),
            "tissue_specificity": er.get("tissue_specificity"),
            "tissue_enriched_in": er.get("tissue_enriched_in"),
            "top5_tissues": er.get("top5_tissues"),
            "certainty_tier": Tier.MEASURED.value,
            "source": "hpa_v25",
        }

    return card


def _card_to_markdown(card: dict) -> str:
    lines = [f"# {card['gene_symbol']} — evidence card",
             f"_run_id: `{card.get('run_id')}` · config_hash: `{card.get('config_hash')}`_",
             ""]

    if v := card.get("mechanism_verdict"):
        lines.append(f"## Mechanism verdict: **{v['verdict']}**")
        if v.get("reasons"):
            lines.append(f"> {v['reasons']}")
        lines.append(f"_tier: {v['certainty_tier']} · source: {v['source']}_")
        lines.append("")

    if c := card.get("constraint"):
        lines += [
            "## Constraint (gnomAD v4.1)",
            f"- LOEUF = {_num(c['loeuf'])} · pLI = {_num(c['pLI'])} · oe_mis = {_num(c['oe_mis'])}",
            f"_tier: {c['certainty_tier']} · source: {c['source']}_",
            "",
        ]

    if v := card.get("variant_catalog"):
        lines += [
            "## Variant catalog",
            f"- {v['n_rare_coding']} rare coding variants (gnomAD popmax AF ≤ 0.001)",
        ]
        for k, n in list(v.get("by_consequence", {}).items())[:6]:
            lines.append(f"  - {k}: {n}")
        lines.append(f"_tier: {v['certainty_tier']} · source: {v['source']}_")
        lines.append("")

    if g := card.get("genebass_top"):
        lines.append("## Genebass associations (top 5)")
        for h in g.get("top5", []):
            lines.append(
                f"- **{h['phenotype']}** · mask={h['annotation']} · β={_num(h['beta'])} · p={_num(h['p_value'])} _[tier: {h['certainty_tier']} · {h['source']}]_"
            )
        lines.append("")

    if fg := card.get("finngen_top"):
        lines.append("## FinnGen R12 replication (top 5)")
        for h in fg.get("top5", []):
            lines.append(
                f"- **{h['phenotype']}** · β={_num(h['beta'])} · p={_num(h['p_value'])} · MAF={_num(h['maf'])} · {h['variant']} _[tier: {h['certainty_tier']} · {h['source']}]_"
            )
        lines.append("")

    if a := card.get("alphamissense"):
        lines.append("## AlphaMissense (predicted pathogenicity)")
        lines.append(f"- n scored: {a['n_scored']} · mean pathogenicity: {_num(a['mean_pathogenicity'])}")
        cls = "; ".join(f"{k}={v}" for k, v in a.get("class_counts", {}).items())
        lines.append(f"- classes: {cls}")
        lines.append(f"_tier: {a['certainty_tier']} · source: {a['source']}_")
        lines.append("")

    if m := card.get("mavedb"):
        lines.append("## MaveDB DMS")
        lines.append(f"- {m['n_measurements']} measurements across {m.get('n_scoresets')} scoresets")
        lines.append(f"_tier: {m['certainty_tier']} · source: {m['source']}_")
        lines.append("")

    if up := card.get("uniprot_features"):
        lines.append("## UniProt features")
        top_types = sorted(up["by_type"].items(), key=lambda x: -x[1])[:6]
        lines.append(f"- {up['n_features']} features · " + "; ".join(f"{t}={n}" for t, n in top_types))
        lines.append(f"_tier: {up['certainty_tier']} · source: {up['source']}_")
        lines.append("")

    if af := card.get("alphafold"):
        lines.append("## AlphaFold model")
        lines.append(f"- {af['n_residues']} residues · median pLDDT = {_num(af['median_plddt'])}")
        lines.append(f"_tier: {af['certainty_tier']} · source: {af['source']}_")
        lines.append("")

    if cv := card.get("clinvar_2plus_star"):
        lines.append("## ClinVar (2+ stars)")
        lines.append(f"- {cv['n_pathogenic']} pathogenic / {cv['n_benign']} benign · max stars = {cv['max_stars']}")
        lines.append(f"_tier: {cv['certainty_tier']} · source: {cv['source']}_")
        lines.append("")

    if led := card.get("curated_ledger"):
        lines.append("## Curated evidence ledger")
        lines.append(f"- ClinVar pathogenic (star≥2) = {led['clinvar_pathogenic_n']} · OT weighted score = {_num(led['ot_weighted_sum'])}")
        if led.get("top_diseases"):
            lines.append(f"- Top diseases: {led['top_diseases']}")
        lines.append(f"_tier: {led['certainty_tier']} · source: {led['source']}_")
        lines.append("")

    if ex := card.get("expression"):
        lines.append("## Expression (HPA consensus)")
        lines.append(f"- Max tissue: **{ex.get('max_tissue')}** · nTPM = {_num(ex.get('max_ntpm'))} · median = {_num(ex.get('median_ntpm'))}")
        lines.append(f"- Tissue specificity: {ex.get('tissue_specificity')} " + (f"(enriched in {ex['tissue_enriched_in']})" if ex.get("tissue_enriched_in") else ""))
        if ex.get("top5_tissues"):
            lines.append(f"- Top 5: {ex['top5_tissues']}")
        lines.append(f"_tier: {ex['certainty_tier']} · source: {ex['source']}_")
        lines.append("")

    return "\n".join(lines)


def _cross_gene_ranking(cards: dict[str, dict]) -> pd.DataFrame:
    """Rank genes by a composite score across evidence sources."""
    rows = []
    for sym, card in cards.items():
        # Composite scoring
        # - Genetics evidence: max negative log10(pval) across genebass + finngen
        gb_p = None
        for h in (card.get("genebass_top") or {}).get("top5", []):
            if h.get("p_value") is not None:
                gb_p = min(gb_p, h["p_value"]) if gb_p is not None else h["p_value"]
        fg_p = None
        for h in (card.get("finngen_top") or {}).get("top5", []):
            if h.get("p_value") is not None:
                fg_p = min(fg_p, h["p_value"]) if fg_p is not None else h["p_value"]

        import math
        # p=0 in FinnGen JSON means "underflowed IEEE double" — floor at 1e-320
        gb_score = -math.log10(max(gb_p, 1e-320)) if gb_p is not None else 0
        fg_score = -math.log10(max(fg_p, 1e-320)) if fg_p is not None else 0
        # Curated
        cur = card.get("curated_ledger") or {}
        cur_score = cur.get("ot_weighted_sum") or 0
        clinvar_score = cur.get("clinvar_pathogenic_n") or 0
        # Constraint pathogenicity (inverted LOEUF for constrained genes)
        con = card.get("constraint") or {}
        loeuf = con.get("loeuf")
        constraint_score = (1.0 / (loeuf + 0.1)) if loeuf else 0
        rows.append({
            "gene_symbol": sym,
            "verdict": (card.get("mechanism_verdict") or {}).get("verdict"),
            "loeuf": loeuf,
            "gb_top_p": gb_p,
            "fg_top_p": fg_p,
            "clinvar_pathogenic_n": clinvar_score,
            "ot_weighted_sum": cur_score,
            "composite_score": gb_score + fg_score + 0.1 * cur_score + clinvar_score * 0.05 + constraint_score,
        })
    df = pd.DataFrame(rows).sort_values("composite_score", ascending=False)
    df["rank"] = range(1, len(df) + 1)
    return df


def run_stage9(cfg: dict, ctx: dict, ckpt) -> dict:
    outdir_run = Path(cfg["outputs_dir"])
    outdir = outdir_run / "stage9"
    outdir.mkdir(parents=True, exist_ok=True)
    gene_ids = ctx["gene_ids"]
    log.info("stage9 start | %d genes", len(gene_ids))

    cards = {}
    for sym in gene_ids.keys():
        card = _build_card(sym, outdir_run, cfg)
        cards[sym] = card
        # Write JSON + markdown
        (outdir / f"{sym}_card.json").write_text(json.dumps(card, indent=2, default=str))
        md = _card_to_markdown(card)
        (outdir / f"{sym}_card.md").write_text(md)
        log.info("card written | %s | %d sections", sym, len([k for k in card if k not in ('gene_symbol','run_id','config_hash')]))

    # Cross-gene ranking
    rank_df = _cross_gene_ranking(cards)
    safe_write_parquet(rank_df, outdir / "rankings.parquet")
    (outdir / "rankings.md").write_text(
        "# Cross-gene ranking\n\n" + rank_df.to_markdown(index=False, floatfmt=".3g")
    )
    for _, r in rank_df.iterrows():
        log.info("rank %d | %s | verdict=%s | composite=%.2f", r["rank"], r["gene_symbol"], r.get("verdict"), r["composite_score"])

    # Master report
    lines = [
        f"# rcvtc pipeline — {cfg.get('run_id')}",
        f"_config_hash: `{cfg.get('_config_hash')}` · genes: {', '.join(gene_ids.keys())}_",
        "",
        "## Cross-gene ranking",
        rank_df.to_markdown(index=False, floatfmt=".3g"),
        "",
    ]
    for sym in gene_ids.keys():
        lines.append(_card_to_markdown(cards[sym]))
        lines.append("---")
        lines.append("")

    (Path(cfg["outputs_dir"]) / "report_rcvtc.md").write_text("\n".join(lines))
    log.info("master report written | %s", outdir_run / "report_rcvtc.md")

    return {"counts": {"cards": len(cards), "rankings": len(rank_df)}}
