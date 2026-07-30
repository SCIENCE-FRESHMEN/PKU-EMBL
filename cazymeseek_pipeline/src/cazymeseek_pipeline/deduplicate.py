"""MMseqs2 可选代表序列簇预处理（本次新增修复5）。

The dbCAN PDF Box 8 gives 95% sequence identity and 95% coverage as its
example for a non-redundant multi-sample protein catalog. This module exposes
those published defaults. Cluster-level abundance is a CAZymeSeek reporting
extension and is labelled as such, not represented as a dbCAN paper result.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def run_mmseqs(mmseqs: str, fasta: Path, output: Path, tmp: Path, identity: float, coverage: float, threads: int) -> Path:
    # 参数直接透传，默认值对应论文 Box 8 的 >0.95 identity 与 >0.95 coverage。
    output.parent.mkdir(parents=True, exist_ok=True); tmp.mkdir(parents=True, exist_ok=True)
    subprocess.run([mmseqs, "easy-cluster", str(fasta), str(output), str(tmp), "--min-seq-id", str(identity), "-c", str(coverage), "--threads", str(threads)], check=True)
    return output.with_name(output.name + "_cluster.tsv")


def read_clusters(cluster_tsv: Path) -> dict[str, str]:
    """读取 MMseqs 两列结果，建立 member -> representative 映射供定量合并使用。"""
    mapping = {}
    with cluster_tsv.open(newline="", encoding="utf-8") as handle:
        for representative, member in csv.reader(handle, delimiter="\t"):
            mapping[member] = representative
    return mapping


def write_representative_fasta(input_fasta: Path, cluster_tsv: Path, output_fasta: Path) -> None:
    """Write nucleotide/protein records whose IDs are MMseqs representatives.

    MMseqs clusters proteins, but the P8 BWA reference is the matching Prokka
    FFN CDS file. Using the same IDs preserves the GFF/FAA/FFN relation required
    by the dbCAN protocol.
    """
    # 仅保留代表 ID 对应的 CDS；ID 不改变，保证 Prokka FAA/FFN/GFF 可追溯。
    representatives = set(read_clusters(cluster_tsv).values())
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    keep = False
    with input_fasta.open(encoding="utf-8") as source, output_fasta.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith(">"):
                keep = line[1:].split()[0] in representatives
            if keep: target.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--fasta", required=True); parser.add_argument("--output", required=True); parser.add_argument("--tmp", required=True); parser.add_argument("--identity", type=float, default=.95); parser.add_argument("--coverage", type=float, default=.95); parser.add_argument("--threads", type=int, default=40); parser.add_argument("--mmseqs", default="mmseqs")
    args = parser.parse_args(); print(run_mmseqs(args.mmseqs, Path(args.fasta), Path(args.output), Path(args.tmp), args.identity, args.coverage, args.threads))
