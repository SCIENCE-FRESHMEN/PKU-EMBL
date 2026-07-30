# CAZymeSeek Annotation Arbitration Refactor Log

## Reference boundary

- Core dbCAN workflow: `2024.01.10.575125v1.full.pdf`, P5--P13 and Box 6--8.
- CAZymeSeek knowledge asset: `EB05_ CAZymeSeek.pdf` and supplied four ChromaDB stores.
- This refactor applies the project's requested `HMMER > Hotpep > DIAMOND` order. dbCAN4's protocol describes dbCAN-sub/eCAMI replacing Hotpep; therefore Hotpep is optional and is only used when a compatible result file is supplied.
- No PUL candidate vote, sequence deep-learning classifier, or fabricated EC assignment is introduced.

## Change matrix

| Defect / required rule | Files | Implementation | Biological rationale | Test |
|---|---|---|---|---|
| Fixed source priority | `annotation.py`, `standardize.py`, `export.py`, `config.example.yaml` | `SOURCE_RANK={hmmer:0, hotpep:1, diamond:2}`; output has `annotation_source` and `source_rank` | Conservatively retain the strongest method where annotations conflict | `test_conflict_priority.py` |
| No blind concatenation | `annotation.py`, `standardize.py` | Individual source parsers feed `resolve_conflicts`; no `pd.concat` or unfiltered union | Prevents lower-priority overlapping calls from being double counted | `test_conflict_priority.py` |
| Multi-domain preservation | `annotation.py`, `standardize.py`, `export.py` | Results are keyed by `gene_id + domain_start + domain_end`; independent intervals all emit | CAZyme fusion proteins can contain multiple functional modules | `test_multidomain.py` |
| >=80% overlap arbitration | `annotation.py`, `config.example.yaml` | `interval_overlap()` uses overlap / shorter interval, configurable `domain_overlap_threshold: 0.80` | Highly overlapping calls represent the same domain; retain lowest source rank | `test_conflict_priority.py` |
| Domain coordinates / family fields | `standardize.py`, `export.py` | Mandatory `domain_start`, `domain_end`, `cazy_family`, `cazy_subfamily` | Separates protein architecture from family/subfamily functional label | `test_multidomain.py` |
| Two-tier EC mapping | `annotation.py`, `standardize.py` | `map_ec()` uses subfamily first, then family, otherwise blank; curated TSV/dbCAN-sub only | Subfamily functional specificity is more precise; no uncurated EC inference | `test_ec_mapping.py` |
| PUL routes unchanged | `pipeline.py`, `standardize.py` | Native `--cgc_substrate` stays; `CGC_substrate_PUL` and `CGC_substrate_vote` are separate read-only exports | dbCAN protocol keeps PUL homology and CAZyme majority vote as two routes | `test_pul_passthrough.py` |
| Optional homolog cluster quantification | `deduplicate.py`, `pipeline.py`, `standardize.py`, config | MMseqs2 representative catalog, configurable 0.95 identity/coverage; default full-CDS mode remains | Protocol Box 8 gives a 95%/95% nonredundant catalog example | `test_deduplication.py` |

## Test commands

```bash
cd cazymeseek_pipeline
PYTHONPATH=src python -m unittest discover -s tests -v
PYTHONPATH=src python demo/test_demo.py
```

The fixture suite checks arbitration, independent domains, EC fallback, unchanged PUL-field pass-through, and representative FASTA construction without requiring external databases.
