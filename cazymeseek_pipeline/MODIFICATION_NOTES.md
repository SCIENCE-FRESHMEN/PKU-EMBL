# CAZymeSeek Pipeline Modification Log

## Scope and evidence boundary

- **Reference:** `2024.01.10.575125v1.full.pdf`, especially P5--P13 and Box 6--8; `EB05_ CAZymeSeek.pdf` plus the supplied four-store vector asset.
- **Boundary:** this record distinguishes published dbCAN behavior from CAZymeSeek reporting extensions. No sequence deep-learning model or new PUL voting algorithm has been introduced.
- **Regression suite:** run `PYTHONPATH=src python -m unittest discover -s tests -v`.

## Issue-to-change matrix

| Teacher issue | Changed files | Implemented logic | Biological/computational basis | Regression test |
|---|---|---|---|---|
| 1. Conflicting dbCAN calls lack HMM priority | `annotation.py`, `standardize.py`, `export.py` | `protein_primary_family` follows `dbCAN HMM > dbCAN-sub/eCAMI > DIAMOND`; raw calls, selected method and conflict flag remain traceable | dbCAN protocol Box 6 gives this family preference | `test_conflict_priority.py` |
| 2. Multi-domain proteins are truncated | `annotation.py`, `standardize.py`, `export.py` | Every HMM/dbCAN-sub hit emits one domain row with index and coordinates | A protein may contain multiple independent CAZyme modules; first-hit-only export loses architecture | `test_multidomain.py` |
| 3. Subfamily EC is missing | `annotation.py`, `standardize.py` | EC only propagates from curated dbCAN-sub/overview fields for the corresponding dbCAN-sub domain | Protocol maps eCAMI subfamilies to EC and substrate; no universal family-suffix EC table is asserted | `test_ec_mapping.py` |
| 4. dbCAN-PUL preference and candidate vote are unchanged | `pipeline.py`, `standardize.py`, `export.py` | Existing `--cgc_substrate` call remains; PUL homology and majority vote export in separate fields; no re-vote code | Protocol prefers dbCAN-PUL homology to majority voting | `test_pul_passthrough.py` |
| 5. Homolog redundancy before optional abundance mapping | `deduplicate.py`, `pipeline.py`, `standardize.py`, `config.example.yaml` | Optional MMseqs representatives; matching FFN is used for BWA, representatives receive one quantitative record and `protein_cluster_id` | Box 8 example uses >0.95 identity and >0.95 coverage. Standard P13 all-CDS `dbcan_utils` path remains default | `test_deduplication.py` |

## Operational notes

- `deduplication.enabled: false` is the default and preserves official P13 aggregate abundance output.
- `deduplication.enabled: true` is a CAZymeSeek cluster-level reporting extension. It must be reported with `protein_cluster_id` and must not be described as the dbCAN paper's mandatory TPM procedure.
- A/B/C fields are a transparent project reporting convention. They are not dbCAN confidence scores.
- Tests use fixtures only. They do not claim synthetic sequences are biologically annotated CAZymes.
