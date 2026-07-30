"""Build domain-resolved standardized annotations from dbCAN tool outputs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .annotation import load_ec_mapping, map_ec, read_tsv, resolve_conflicts, source_candidates, field
from .export import export_rows


def index_rows(path: Path | None, key: str) -> dict[str, dict[str, str]]:
    """建立外部表索引；不对不同来源的注释表做 pd.concat 式盲合并。"""
    if not path or not path.exists():
        return {}
    return {value: row for row in read_tsv(path) if (value := field(row, key))}


def standardize(sample_id: str, abundance: Path, output: Path, hmmer: Path | None = None,
                hotpep: Path | None = None, diamond: Path | None = None,
                dbcan_sub: Path | None = None, ec_mapping: Path | None = None,
                cgc_standard: Path | None = None, substrate: Path | None = None,
                cluster_map: Path | None = None, domain_overlap_threshold: float = .80) -> None:
    """输出逐域、可审计结果。

    规则：HMMER > Hotpep > DIAMOND，按 gene_id + 域区间判重；同一蛋白
    的独立域均保留。dbCAN-sub仅提供策展 EC/底物映射。PUL 两路线只透传。
    """
    raw = source_candidates(hmmer, "hmmer") + source_candidates(hotpep, "hotpep") + source_candidates(diamond, "diamond")
    accepted = resolve_conflicts(raw, domain_overlap_threshold)
    abundance_by_gene = index_rows(abundance, "protein_sequence_id")
    cgc_by_gene = index_rows(cgc_standard, "gene_id")
    substrate_by_cgc = index_rows(substrate, "cgc_id")
    # dbCAN-sub HMM output can carry EC per hit. 自定义映射表优先；没有表时读取该策展输出。
    ec_by_subfamily = load_ec_mapping(ec_mapping or dbcan_sub)
    clusters = {}
    if cluster_map and cluster_map.exists():
        with cluster_map.open(newline="", encoding="utf-8") as handle:
            clusters = {member: representative for representative, member in csv.reader(handle, delimiter="\t")}
    output_rows = []
    emitted_clusters: set[tuple[str, str, int | None, int | None]] = set()
    for candidate in accepted:
        gene = str(candidate["gene_id"]); representative = clusters.get(gene, gene)
        # 去冗余模式：同一代表序列的同一区间仅计数一次；不同域仍完整保留。
        cluster_key = (representative, str(candidate["cazy_subfamily"]), candidate["domain_start"], candidate["domain_end"])
        if clusters and cluster_key in emitted_clusters:
            continue
        emitted_clusters.add(cluster_key)
        counts = abundance_by_gene.get(representative, abundance_by_gene.get(gene, {}))
        cgc = cgc_by_gene.get(gene, {}); cgc_id = field(cgc, "cgc_id"); substrate_row = substrate_by_cgc.get(cgc_id, {})
        output_rows.append({
            "sample_id": sample_id, "contig_id": field(cgc, "contig_id"), "gene_id": gene,
            "protein_sequence_id": gene, "protein_cluster_id": representative,
            "domain_start": candidate["domain_start"] or "", "domain_end": candidate["domain_end"] or "",
            "cazy_family": candidate["cazy_family"], "cazy_subfamily": candidate["cazy_subfamily"],
            "annotation_source": candidate["annotation_source"], "source_rank": candidate["source_rank"],
            "EC": map_ec(candidate, ec_by_subfamily), "annotation_evidence": candidate["raw_annotation"],
            "CGC_id": cgc_id, "CGC_gene_composition": field(cgc, "annotation"),
            # 问题8：原始 dbCAN-PUL 同源与多数投票结果只读、分列保存，不重投票。
            "CGC_substrate_PUL": field(substrate_row, "substrate_of_the_hit_pul"),
            "CGC_substrate_vote": field(substrate_row, "substrate_predicted_by_majority_voting_of_cazymes_in_cgc"),
            "TPM": counts.get("TPM", ""), "RPM": counts.get("RPM", ""), "RPKM": counts.get("RPKM", ""),
        })
    export_rows(output_rows, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--sample-id", required=True); parser.add_argument("--abundance", required=True); parser.add_argument("--output", required=True); parser.add_argument("--hmmer"); parser.add_argument("--hotpep"); parser.add_argument("--diamond"); parser.add_argument("--dbcan-sub"); parser.add_argument("--ec-mapping"); parser.add_argument("--cgc-standard"); parser.add_argument("--substrate"); parser.add_argument("--cluster-map"); parser.add_argument("--domain-overlap-threshold", type=float, default=.80)
    args = parser.parse_args()
    standardize(args.sample_id, Path(args.abundance), Path(args.output), *(Path(value) if value else None for value in (args.hmmer, args.hotpep, args.diamond, args.dbcan_sub, args.ec_mapping, args.cgc_standard, args.substrate, args.cluster_map)), args.domain_overlap_threshold)
