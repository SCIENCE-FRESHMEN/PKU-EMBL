"""Published dbCAN read-to-annotation workflow wrapper.

Source convention: ``[dbCAN PDF, section]`` refers to
2024.01.10.575125v1.full.pdf. This module only orchestrates published tools;
run_dbcan performs family, subfamily, CGC and substrate annotation.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import yaml


def run(command: list[str], log: Path) -> None:
    """Execute one auditable command and capture stdout/stderr in its stage log."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)


def sample_rows(path: Path):
    """Read ``sample_id, read1, read2`` paired-end input rows."""
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def trim_output(read: str, mate: int) -> str:
    """Return Trim Galore paired-end output name [dbCAN PDF, P2]."""
    name = Path(read).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)] + f"_val_{mate}.fq.gz"
    raise ValueError(f"Unsupported FASTQ suffix: {read}")


def dbcan_command(cfg: dict, faa: Path, gff: Path, output: Path) -> list[str]:
    """Use P7 syntax: all three CAZyme methods, CGC finding and two substrate routes.

    No ``--tools`` option is passed: the protocol states this activates dbCAN
    HMMER, DIAMOND/CAZy and dbCAN-sub HMMER. ``--cgc_substrate`` causes both
    dbCAN-PUL homology and dbCAN-sub majority-voting predictions [P5--P7].
    """
    cpu = str(cfg["resources"]["threads"])
    return [cfg["tools"]["run_dbcan"], str(faa), "protein", "--db_dir", cfg["paths"]["dbcan_db_dir"],
            "--out_dir", str(output), "--dia_cpu", cpu, "--hmm_cpu", cpu,
            "--dbcan_thread", cpu, "--tf_cpu", cpu, "--stp_cpu", cpu,
            "-c", str(gff), "--cgc_substrate"]


def individual_commands(cfg: dict, sample: dict, root: Path) -> list[tuple[str, list[str]]]:
    """P1--P7 individual-assembly commands [dbCAN PDF, Fig. 3; P1--P7]."""
    sid, r1, r2 = sample["sample_id"], sample["read1"], sample["read2"]
    out = root / "output" / sid
    clean, assembly, prokka, dbcan = (out / part for part in ("clean", "assembly", "prokka", "dbcan"))
    t1, t2 = clean / trim_output(r1, 1), clean / trim_output(r2, 2)
    contigs = assembly / f"{sid}.contigs.fa"
    threads = str(cfg["resources"]["threads"])
    commands = [
        # P1 assesses possible contamination; reference alignment removal is optional below.
        ("kraken", [cfg["tools"]["kraken2"], "--threads", threads, "--quick", "--paired", "--db", cfg["paths"]["kraken_db_dir"], "--report", str(clean / f"{sid}.kreport"), "--output", str(clean / f"{sid}.kraken.out"), r1, r2]),
        # P2 retains only paired reads; no unpaired-read option is used in the protocol.
        ("trim", [cfg["tools"]["trim_galore"], "--paired", "--illumina", "--cores", threads, "--output_dir", str(clean), r1, r2]),
    ]
    reads = (t1, t2)
    removal = cfg.get("contamination_removal", {})
    if removal.get("enabled", False):
        # P1 Box 1 removes pairs unmapped to an explicitly selected contamination reference.
        prefix = clean / "contaminant"
        unmapped = clean / "unmapped.bam"
        reads = (clean / "clean_1.fq.gz", clean / "clean_2.fq.gz")
        commands.extend([
            ("contaminant_index", [cfg["tools"]["bwa"], "index", "-p", str(prefix), removal["reference_fasta"]]),
            ("contaminant_map", [cfg["tools"]["bwa"], "mem", "-t", threads, str(prefix), str(t1), str(t2)]),
        ])
        # This special pipeline is executed in ``process_individual`` to preserve BWA -> SAMtools streaming.
    commands.extend([
        # P3 uses 1,000 bp minimum contigs and 50% available memory in the published example.
        ("assembly", [cfg["tools"]["megahit"], "-m", str(cfg["assembly"]["megahit_memory_fraction"]), "-t", threads, "-o", str(assembly), "--out-prefix", sid, "--min-contig-len", str(cfg["assembly"]["min_contig_len"]), "-1", str(reads[0]), "-2", str(reads[1])]),
        # P4: bacterial Prokka gene calls provide matched FAA, FFN and GFF identifiers.
        ("prokka", [cfg["tools"]["prokka"], "--kingdom", "Bacteria", "--cpus", threads, "--outdir", str(prokka), "--prefix", sid, "--addgenes", "--addmrna", "--locustag", sid, str(contigs)]),
        ("dbcan", dbcan_command(cfg, prokka / f"{sid}.faa", prokka / f"{sid}.gff", dbcan)),
    ])
    return commands


def process_individual(cfg: dict, sample: dict) -> None:
    """Run P1--P13, including CDS and contig mapping plus dbcan_utils abundance."""
    root = Path(cfg["paths"]["project_dir"]); sid = sample["sample_id"]
    out = root / "output" / sid; logs = root / "output" / "logs" / sid
    for name, command in individual_commands(cfg, sample, root):
        if name == "contaminant_map":
            # Deliberately fail closed: optional removal needs its paired-output stream implemented for the local SAMtools release.
            raise RuntimeError("contamination_removal.enabled requires a locally validated paired-end BWA/SAMtools extraction command; see README Box 1.")
        run(command, logs / f"{name}.log")
    prokka, assembly, abundance, dbcan = out / "prokka", out / "assembly", out / "abundance", out / "dbcan"
    abundance.mkdir(parents=True, exist_ok=True)
    clean = out / "clean"; t1, t2 = clean / trim_output(sample["read1"], 1), clean / trim_output(sample["read2"], 2)
    threads = str(cfg["resources"]["threads"])
    # P8--P11: map clean reads to all CDS, count with Bedtools, then normalize with dbcan_utils.
    for label, reference, bam_name in [("CDS", prokka / f"{sid}.ffn", f"{sid}.CDS.bam"), ("contig", assembly / f"{sid}.contigs.fa", f"{sid}.contig.bam")]:
        run([cfg["tools"]["bwa"], "index", str(reference)], logs / f"{label}.index.log")
        sam = out / f"{sid}.{label}.sam"; bam = out / bam_name
        run([cfg["tools"]["bwa"], "mem", "-t", threads, "-o", str(sam), str(reference), str(t1), str(t2)], logs / f"{label}.map.log")
        run([cfg["tools"]["samtools"], "sort", "-@", threads, "-o", str(bam), str(sam)], logs / f"{label}.sort.log")
    # P11 needs a two-column genome-length file and a BED-like zero-start record per CDS.
    lengths = subprocess.check_output([cfg["tools"]["seqkit"], "fx2tab", "-l", "-n", "-i", str(prokka / f"{sid}.ffn")], text=True)
    length_file, bed_file = abundance / f"{sid}.length", abundance / f"{sid}.bed"
    with length_file.open("w", encoding="utf-8") as length_handle, bed_file.open("w", encoding="utf-8") as bed_handle:
        for line in lengths.splitlines():
            gene, length = line.split("\t")[:2]
            length_handle.write(f"{gene}\t{length}\n")
            bed_handle.write(f"{gene}\t0\t{length}\n")
    # Bedtools writes coverage to stdout; capture it as the P11 depth input.
    with (abundance / f"{sid}.depth.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([cfg["tools"]["bedtools"], "coverage", "-g", str(length_file), "-sorted", "-a", str(bed_file), "-counts", "-b", str(out / f"{sid}.CDS.bam")], check=True, stdout=handle, stderr=subprocess.STDOUT)
    # dbcan_utils is the published normalizer; it computes gene/family/subfamily/CGC/substrate TPM.
    for utility in ("fam_abund", "fam_substrate_abund", "CGC_abund", "CGC_substrate_abund"):
        run([cfg["tools"]["dbcan_utils"], utility, "-bt", str(abundance / f"{sid}.depth.txt"), "-i", str(dbcan), "-a", "TPM"], logs / f"{utility}.log")
    # The fixed project schema is a sidecar; published dbcan_utils aggregate files remain authoritative.
    run(["python", "-m", "cazymeseek_pipeline.abundance", "--bam", str(out / f"{sid}.CDS.bam"), "--ffn", str(prokka / f"{sid}.ffn"), "--output", str(abundance / "gene_abundance.tsv")], logs / "gene_abundance.log")
    run(["python", "-m", "cazymeseek_pipeline.standardize", "--sample-id", sid, "--overview", str(dbcan / "overview.txt"), "--abundance", str(abundance / "gene_abundance.tsv"), "--cgc-standard", str(dbcan / "cgc_standard.out"), "--substrate", str(dbcan / "substrate.out"), "--output", str(out / "standardized_annotations.csv")], logs / "standardize.log")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True); parser.add_argument("--sample")
    args = parser.parse_args()
    with open(args.config, encoding="utf-8") as handle: cfg = yaml.safe_load(handle)
    mode = cfg["assembly"]["mode"]
    if mode != "individual":
        raise NotImplementedError(f"'{mode}' is documented in README but not silently approximated by this individual-assembly wrapper.")
    for sample in sample_rows(Path(cfg["paths"]["samples_tsv"])):
        if not args.sample or sample["sample_id"] == args.sample: process_individual(cfg, sample)


if __name__ == "__main__": main()
