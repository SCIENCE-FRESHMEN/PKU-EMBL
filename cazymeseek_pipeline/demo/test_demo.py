"""Offline minimal demo: validates schema/export and renders protocol-style figures.

It intentionally does not claim that synthetic reads yield biological CAZyme calls.
The full read-to-annotation path requires Linux bioinformatics dependencies and
the downloaded dbCAN database described in the companion protocol.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cazymeseek_pipeline.export import export_rows
from cazymeseek_pipeline.visualize import barplot, heatmap


def main() -> None:
    output = ROOT / "output" / "demo"; output.mkdir(parents=True, exist_ok=True)
    # Small input artefacts demonstrate accepted FASTQ/FAA forms without annotation claims.
    (output / "demo_R1.fastq").write_text("@demo/1\nACGTACGT\n+\nFFFFFFFF\n", encoding="ascii")
    (output / "demo_R2.fastq").write_text("@demo/2\nTGCATGCA\n+\nFFFFFFFF\n", encoding="ascii")
    (output / "demo.faa").write_text(">demo_gene\nMKKLLVL\n", encoding="ascii")
    # Synthetic output verifies all fixed columns and the transparent A/B/C export rules.
    export_rows([
        {"sample_id": "demo", "gene_id": "gene_A", "protein_sequence_id": "gene_A", "family_id": "GH1", "subfamily_id": "GH1_1", "EC": "3.2.1.21", "annotation_evidence": "dbCAN; dbCAN_sub", "TPM": 100},
        {"sample_id": "demo", "gene_id": "gene_B", "protein_sequence_id": "gene_B", "family_id": "GH43", "annotation_evidence": "dbCAN", "TPM": 40},
        {"sample_id": "demo", "gene_id": "gene_C", "protein_sequence_id": "gene_C", "CGC_substrate_vote": "xylan", "annotation_evidence": "majority_vote", "TPM": 5},
    ], output / "standardized_annotations.csv")
    with (output / "fam_abund.out").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["Family", "TPM"]); writer.writerows([["GH1", 100], ["GH43", 40]])
    with (output / "fam_substrate_abund.out").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["Substrate", "TPM"]); writer.writerows([["cellobiose", 100], ["xylan", 40]])
    barplot(str(output / "fam_abund.out"), str(output / "family_abundance"), "Family", "TPM")
    heatmap(str(output / "fam_substrate_abund.out"), str(output / "substrate_heatmap"))
    print(f"Demo completed: {output}")


if __name__ == "__main__": main()
