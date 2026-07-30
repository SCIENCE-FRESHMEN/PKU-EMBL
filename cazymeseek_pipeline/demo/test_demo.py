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
from cazymeseek_pipeline.standardize import standardize
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
    # NEW fixes 1--3: conflict resolves to dbCAN HMM, two independent domains
    # persist as two rows, and EC is read from dbCAN-sub/overview data.
    (output / "overview.txt").write_text("Gene ID\tEC#\tdbCAN\tdbCAN_sub\tDIAMOND\nprotein_1\t3.2.1.21\tGH1(5-90);GH3(110-200)\tGH43_e149\tGH43_4\n", encoding="utf-8")
    (output / "hmmer.out").write_text("Gene ID\tHMM Profile\tAli From\tAli To\nprotein_1\tGH1\t5\t90\nprotein_1\tGH3\t110\t200\n", encoding="utf-8")
    (output / "dbcan-sub.hmm.out").write_text("Gene ID\tHMM Profile\tEC\tSubstrate\nprotein_1\tGH43_e149\t3.2.1.21\txylan\n", encoding="utf-8")
    (output / "gene_abundance.tsv").write_text("protein_sequence_id\tTPM\tRPM\tRPKM\nprotein_1\t100\t50\t50\n", encoding="utf-8")
    standardize("demo", output / "overview.txt", output / "gene_abundance.tsv", output / "domain_resolved.csv", output / "hmmer.out", output / "dbcan-sub.hmm.out")
    with (output / "domain_resolved.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3, "Multiple HMM/dbCAN-sub domains must not be truncated."
    assert rows[0]["protein_primary_family"] == "GH1" and rows[0]["selected_method"] == "dbCAN"
    assert {row["family_id"] for row in rows} == {"GH1", "GH3", "GH43"}
    assert any(row["EC"] == "3.2.1.21" for row in rows), "EC must propagate from dbCAN-sub/overview."
    print(f"Demo completed: {output}")


if __name__ == "__main__": main()
