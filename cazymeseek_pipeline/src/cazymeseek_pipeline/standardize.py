"""Adapt documented run_dbcan overview.txt and dbcan_utils TPM outputs to one schema."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from .export import export_rows

FAMILY = re.compile(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+(?:_\d+)?)\b")


def column(row: dict[str, str], *names: str) -> str:
    lowered = {key.lower().replace(" ", "_"): value for key, value in row.items()}
    for name in names:
        if name in lowered:
            return lowered[name]
    return ""


def family_and_subfamily(annotation: str) -> tuple[str, str]:
    tokens = FAMILY.findall(annotation or "")
    if not tokens:
        return "", ""
    first = tokens[0]
    return first.split("_")[0], first if "_" in first else ""


def standardize(sample_id: str, overview: Path, abundance: Path, output: Path, cgc_standard: Path | None = None, substrate: Path | None = None) -> None:
    """Parse overview.txt fields documented in Box 6 of the dbCAN PDF.

    The PDF says dbCAN > dbCAN-sub/eCAMI > DIAMOND when families disagree;
    this adapter therefore selects the first dbCAN-family token before later
    tokens in the ordered overview row. It does not invent a consensus call.
    """
    with abundance.open(newline="", encoding="utf-8") as handle:
        tpm = {column(row, "protein_sequence_id", "gene_id"): row for row in csv.DictReader(handle, delimiter="\t")}
    cgc_by_gene: dict[str, dict[str, str]] = {}
    if cgc_standard and cgc_standard.exists():
        with cgc_standard.open(newline="", encoding="utf-8") as handle:
            for source in csv.DictReader(handle, delimiter="\t"):
                gene = column(source, "gene_id")
                if gene: cgc_by_gene[gene] = source
    substrate_by_cgc: dict[str, dict[str, str]] = {}
    if substrate and substrate.exists():
        with substrate.open(newline="", encoding="utf-8") as handle:
            for source in csv.DictReader(handle, delimiter="\t"):
                cgc = column(source, "cgc_id")
                if cgc: substrate_by_cgc[cgc] = source
    rows = []
    with overview.open(newline="", encoding="utf-8") as handle:
        for source in csv.DictReader(handle, delimiter="\t"):
            gene = column(source, "gene_id", "gene_id#", "protein_id")
            evidence = " | ".join(value for value in source.values() if value and value != "-")
            family, subfamily = family_and_subfamily(evidence)
            counts = tpm.get(gene, {})
            cgc = cgc_by_gene.get(gene, {}); cgc_id = column(cgc, "cgc_id")
            substrate_row = substrate_by_cgc.get(cgc_id, {})
            rows.append({"sample_id": sample_id, "gene_id": gene, "protein_sequence_id": gene,
                         "contig_id": column(cgc, "contig_id"), "cazy_class": family[:2],
                         "family_id": family, "subfamily_id": subfamily,
                         "EC": column(source, "ec#", "ec"), "CGC_id": cgc_id,
                         "CGC_gene_composition": column(cgc, "annotation"),
                         "CAZyme_substrate": column(source, "substrate"),
                         "CGC_substrate_PUL": column(substrate_row, "substrate_of_the_hit_pul", "substrate_of_the_hit_pul_"),
                         "CGC_substrate_vote": column(substrate_row, "substrate_predicted_by_majority_voting_of_cazymes_in_cgc"),
                         "annotation_evidence": evidence, "TPM": counts.get("TPM", ""),
                         "RPM": counts.get("RPM", ""), "RPKM": counts.get("RPKM", "")})
    export_rows(rows, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--sample-id", required=True); parser.add_argument("--overview", required=True); parser.add_argument("--abundance", required=True); parser.add_argument("--output", required=True); parser.add_argument("--cgc-standard"); parser.add_argument("--substrate")
    args = parser.parse_args(); standardize(args.sample_id, Path(args.overview), Path(args.abundance), Path(args.output), Path(args.cgc_standard) if args.cgc_standard else None, Path(args.substrate) if args.substrate else None)
