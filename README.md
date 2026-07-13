# RC-VTC: Rare-Coding-Variant Target Characterization

Reusable Python pipeline for characterizing drug target–disease relationships in a small (<20 gene) target list, using public rare-coding-variant human genetics.

Full design and rationale: `execution_trace/PLAN.md`.

## Install

```bash
uv pip install -e .
```

## Run

```bash
python run_pipeline.py --config config/example_config.yaml
```

Outputs land in `/mnt/results/rcvtc_runs/<run_id>/`.

## Certainty tiers

Every claim in the output is tagged one of:

- **measured** — experimental data (MaveDB DMS, ClinVar PS3/BS3, gnomAD population observations)
- **inferred** — statistical test on summary statistics (Genebass, ExPheWas, OT, FinnGen)
- **predicted** — computational model (AlphaMissense, ESM, REVEL, LOEUF-anchored predictions)
