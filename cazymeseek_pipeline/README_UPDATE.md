# CAZymeSeek Pipeline: Five-Item Update Guide

## What changed

| Update | User-visible behavior | Output fields |
|---|---|---|
| Method conflicts | Protein primary family follows dbCAN HMM, then dbCAN-sub/eCAMI, then DIAMOND | `protein_primary_family`, `selected_method`, `annotation_conflict`, `annotation_evidence` |
| Multi-domain proteins | Each domain is a separate row; no first-domain truncation | `domain_index`, `domain_start`, `domain_end`, domain-level `family_id`, `subfamily_id` |
| EC mapping | EC is retained only when supplied by curated dbCAN-sub/overview evidence | `EC`, `CAZyme_substrate` |
| PUL routes | dbCAN-PUL homology and majority-vote calls stay separate and unchanged | `CGC_substrate_PUL`, `CGC_substrate_vote` |
| Optional homolog reduction | MMseqs representative CDS catalog can be used before BWA quantification | `protein_cluster_id`, TPM/RPM/RPKM |

## Standard run

```bash
cd /data/cazymeseek_pipeline
conda activate cazymeseek
export PYTHONPATH="$PWD/src"
python -m cazymeseek_pipeline.pipeline --config config/config.yaml --sample SAMPLE_ID
```

`run_dbcan --cgc_substrate` is retained in the generated command. It produces the original `substrate.out`; the pipeline only exports its two reported methods without altering their candidate vote.

## Optional MMseqs representative-cluster quantification

In `config/config.yaml`:

```yaml
deduplication:
  enabled: true
  min_seq_identity: 0.95
  coverage: 0.95
```

- Keep `enabled: false` to obtain the original full-CDS `dbcan_utils` P13 family/CGC/substrate abundance path.
- Use `enabled: true` when the study requires representative-cluster quantification to reduce repeated homolog measurement. The output contains `protein_cluster_id`; interpret TPM at the representative cluster level.
- The defaults match the protocol Box 8 non-redundant protein catalog example. Do not alter the values without recording the study-specific rationale.

## Interpretation of the updated fields

- `protein_primary_family`: resolved protein-level family only; conflict priority is applied here.
- `family_id` and `subfamily_id`: domain-level labels. A fusion protein can have more than one result row.
- `EC`: only a curated dbCAN-sub/overview EC. Empty means no curated mapping was found, not that an EC was inferred as absent.
- `annotation_evidence`: all raw dbCAN/dbCAN-sub/DIAMOND calls retained for review.
- `CGC_substrate_PUL`: dbCAN-PUL homology result, preferred by the published protocol when present.
- `CGC_substrate_vote`: dbCAN-sub CAZyme majority-voting result; it remains a separate route.

## Automated regression tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The suite independently checks all five changes with local fixtures; it does not need run_dbcan databases or Linux bioinformatics executables.
