"""Normalize run_dbcan output and attach conservative evidence tiers."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

STANDARD_FIELDS = ["sample_id", "contig_id", "gene_id", "protein_sequence_id", "cazy_class", "family_id", "subfamily_id", "EC", "annotation_evidence", "CGC_id", "CGC_gene_composition", "CAZyme_substrate", "CGC_substrate_PUL", "CGC_substrate_vote", "TPM", "RPM", "RPKM", "confidence_tier", "source_database", "version"]


def tier(row: dict[str, str]) -> str:
    """Project reporting convention, not a dbCAN confidence score.

    The PDF gives a preference order (dbCAN-PUL homology > majority voting) and
    a family-method preference (dbCAN > dbCAN-sub > DIAMOND), but it does not
    define A/B/C labels. This transparent mapping preserves those evidence types.
    """
    if row.get("CGC_substrate_PUL") or (row.get("subfamily_id") and row.get("EC")):
        return "A"
    if row.get("family_id") and row.get("annotation_evidence"):
        return "B"
    return "C"


def export_rows(rows: list[dict[str, Any]], output: str | Path, version: str = "run_dbcan-4.0.0") -> None:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STANDARD_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for source in rows:
            row = {field: str(source.get(field, "")) for field in STANDARD_FIELDS}
            row["confidence_tier"] = tier(row)
            row["source_database"] = row["source_database"] or "dbCAN/CAZy/dbCAN-PUL"
            row["version"] = row["version"] or version
            writer.writerow(row)
