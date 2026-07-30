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
        {"sample_id": "demo", "gene_id": "gene_A", "protein_sequence_id": "gene_A", "cazy_family": "GH1", "cazy_subfamily": "GH1_1", "EC": "3.2.1.21", "annotation_source": "hmmer", "source_rank": 0, "annotation_evidence": "HMMER:GH1", "TPM": 100},
        {"sample_id": "demo", "gene_id": "gene_B", "protein_sequence_id": "gene_B", "cazy_family": "GH43", "annotation_source": "diamond", "source_rank": 2, "annotation_evidence": "DIAMOND:GH43", "TPM": 40},
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
    (output / "hmmer.out").write_text("Gene ID\tHMM Profile\tAli From\tAli To\nprotein_1\tGH1\t5\t90\nprotein_1\tGH3\t110\t200\n", encoding="utf-8")
    (output / "hotpep.out").write_text("Gene ID\tHMM Profile\tAli From\tAli To\nprotein_1\tGH43_4\t6\t89\n", encoding="utf-8")
    (output / "diamond.out").write_text("Gene ID\tHMM Profile\tAli From\tAli To\nprotein_1\tGH43_4\t6\t89\n", encoding="utf-8")
    (output / "ec_mapping.tsv").write_text("subfamily\tEC\nGH43_4\t3.2.1.21\nGH43\t3.2.1.99\n", encoding="utf-8")
    (output / "gene_abundance.tsv").write_text("protein_sequence_id\tTPM\tRPM\tRPKM\nprotein_1\t100\t50\t50\n", encoding="utf-8")
    standardize("demo", output / "gene_abundance.tsv", output / "domain_resolved.csv", hmmer=output / "hmmer.out", hotpep=output / "hotpep.out", diamond=output / "diamond.out", ec_mapping=output / "ec_mapping.tsv")
    with (output / "domain_resolved.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2, "HMMER 的两个不重叠结构域必须保留，重叠 Hotpep/DIAMOND 必须裁决掉。"
    assert {row["cazy_family"] for row in rows} == {"GH1", "GH3"}
    assert {row["source_rank"] for row in rows} == {"0"}
    print(f"Demo completed: {output}")


if __name__ == "__main__": main()
