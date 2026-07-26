"""CDS-level RPM, RPKM and TPM sidecar calculated from BWA/SAMtools summaries.

The protocol's authoritative aggregate values are calculated by dbcan_utils from
Bedtools P11 counts. This sidecar carries the same documented normalization
formula into the standardized per-gene export [dbCAN PDF, P11--P13].
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths, current, size = {}, None, 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current: lengths[current] = size
            current, size = line[1:].split()[0], 0
        else: size += len(line.strip())
    if current: lengths[current] = size
    return lengths


def quantify(bam: Path, ffn: Path, output: Path) -> None:
    """Use mapped primary alignments per CDS; multi-mapping policy must be reported downstream."""
    lengths = fasta_lengths(ffn)
    idxstats = subprocess.check_output(["samtools", "idxstats", str(bam)], text=True)
    records = []
    for line in idxstats.splitlines():
        reference, _, mapped, _ = line.split("\t")
        if reference in lengths and lengths[reference]:
            records.append((reference, lengths[reference], int(mapped)))
    library = sum(record[2] for record in records) or 1
    rpk = {gene: mapped / (length / 1000) for gene, length, mapped in records}
    scale = sum(rpk.values()) or 1
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_sequence_id", "length_bp", "mapped_reads", "RPM", "RPKM", "TPM"], delimiter="\t")
        writer.writeheader()
        for gene, length, mapped in records:
            writer.writerow({"protein_sequence_id": gene, "length_bp": length, "mapped_reads": mapped,
                             "RPM": mapped * 1_000_000 / library,
                             "RPKM": mapped * 1_000_000_000 / (length * library),
                             "TPM": rpk[gene] * 1_000_000 / scale})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--bam", required=True); parser.add_argument("--ffn", required=True); parser.add_argument("--output", required=True)
    args = parser.parse_args(); quantify(Path(args.bam), Path(args.ffn), Path(args.output))
