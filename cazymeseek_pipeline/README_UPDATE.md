# CAZymeSeek Pipeline: Domain-Arbitration Update Guide

## What changed

| Update | User-visible behavior | Output fields |
|---|---|---|
| Method conflicts | Same-domain calls follow HMMER, then Hotpep, then DIAMOND | `annotation_source`, `source_rank`, `annotation_evidence` |
| Multi-domain proteins | Each independent domain is a separate row; >=80% overlap is arbitrated | `domain_start`, `domain_end`, `cazy_family`, `cazy_subfamily` |
| EC mapping | Curated mapping uses subfamily first, then family; unmatched EC stays blank | `EC` |
| PUL routes | dbCAN-PUL homology and majority-vote calls stay separate and unchanged | `CGC_substrate_PUL`, `CGC_substrate_vote` |
| Optional homolog reduction | MMseqs representative CDS catalog can be used before BWA quantification | `protein_cluster_id`, TPM/RPM/RPKM |

## Standard run

```bash
cd /data/cazymeseek_pipeline
conda activate cazymeseek
export PYTHONPATH="$PWD/src"
python -m cazymeseek_pipeline.pipeline --config config/config.yaml --sample SAMPLE_ID
```

`run_dbcan --cgc_substrate` is retained in the generated command. It produces the original `substrate.out`; the pipeline only exports its two reported methods without altering their candidate vote. `annotation.hotpep_output` is optional because modern dbCAN's documented workflow uses dbCAN-sub/eCAMI; when a legacy Hotpep table is supplied it participates with rank 1.

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

- `cazy_family` / `cazy_subfamily`: domain-level family labels. A fusion protein can have more than one result row.
- `domain_start` / `domain_end`: protein-coordinate interval for each output domain.
- `annotation_source` / `source_rank`: HMMER/0, Hotpep/1, DIAMOND/2.
- `EC`: curated subfamily mapping first, family fallback second. Empty means no curated mapping was found.
- `annotation_evidence`: source annotation retained for review.
- `CGC_substrate_PUL`: dbCAN-PUL homology result, preferred by the published protocol when present.
- `CGC_substrate_vote`: dbCAN-sub CAZyme majority-voting result; it remains a separate route.

## Automated regression tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The suite independently checks all five changes with local fixtures; it does not need run_dbcan databases or Linux bioinformatics executables.
