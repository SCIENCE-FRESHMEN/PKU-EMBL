"""从 run_dbcan 原始表生成可审计、逐域展开的 CAZymeSeek 结果。"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from .annotation import domain_rows, ec_values, field, preferred_family, read_tsv
from .export import export_rows


def index_rows(path: Path | None, key: str) -> dict[str, dict[str, str]]:
    """按指定字段建立索引；缺失的可选 dbCAN 文件保持为空索引。"""
    if not path or not path.exists(): return {}
    result = {}
    for row in read_tsv(path):
        value = field(row, key)
        if value: result[value] = row
    return result


def standardize(sample_id: str, overview: Path, abundance: Path, output: Path,
                hmmer: Path | None = None, dbcan_sub: Path | None = None,
                cgc_standard: Path | None = None, substrate: Path | None = None,
                cluster_map: Path | None = None) -> None:
    """Export one row per recognized domain, retaining conflict/audit fields.

    NEW fix 1: the selected family uses dbCAN HMM > dbCAN-sub > DIAMOND.
    NEW fix 2: all domain rows from HMM/domain outputs are retained.
    NEW fix 3: EC values only come from dbCAN-sub/overview curated outputs.
    FIX 4: substrate.out is passed through without changing PUL voting fields.
    FIX 5: optional MMseqs representative IDs attach to quantitative records.
    """
    # 【原有逻辑】读取丰度、CGC 与底物表；此处不重算 PUL 投票。
    abundance_by_gene = index_rows(abundance, "protein_sequence_id")
    cgc_by_gene = index_rows(cgc_standard, "gene_id")
    substrate_by_cgc = index_rows(substrate, "cgc_id")
    clusters = {}
    # 【新增修复5】开启去冗余时，用 member -> representative 映射读取代表簇 TPM。
    if cluster_map and cluster_map.exists():
        with cluster_map.open(newline="", encoding="utf-8") as handle:
            clusters = {member: representative for representative, member in csv.reader(handle, delimiter="\t")}
    domains_by_gene = domain_rows(hmmer, dbcan_sub)
    output_rows = []
    for overview_row in read_tsv(overview):
        gene = field(overview_row, "gene_id", "gene_id_", "protein_id")
        if not gene: continue
        # 【新增修复5】每个代表簇只输出一次定量记录，避免同一代表序列的 reads 重复写入成员基因。
        if clusters and gene != clusters.get(gene, gene):
            continue
        family, selected_method, calls = preferred_family(overview_row)
        # 【兼容原有逻辑】若 HMM 明细文件不存在，保留 overview 层一行，避免有效调用丢失。
        domains = domains_by_gene.get(gene, [{"method": selected_method, "domain_annotation": family,
                                              "domain_start": "", "domain_end": "", "domain_ec": "", "domain_substrate": ""}])
        counts = abundance_by_gene.get(clusters.get(gene, gene), abundance_by_gene.get(gene, {}))
        cgc = cgc_by_gene.get(gene, {}); cgc_id = field(cgc, "cgc_id"); substrate_row = substrate_by_cgc.get(cgc_id, {})
        evidence = " | ".join(f"{method}:{','.join(values) or '-'}" for method, values in calls.items())
        conflict = len({value.split("_")[0] for values in calls.values() for value in values}) > 1
        for index, domain in enumerate(domains, start=1):
            domain_families = domain["domain_annotation"].split()
            domain_family = next((value.split("_")[0] for value in domain_families if value.startswith(("GH", "GT", "PL", "CE", "AA", "CBM"))), family)
            # 【新增修复1+2】protein_primary_family 遵循 HMM 优先；但 domain family_id
            # 保持独立，避免把 GH1--GH3 融合蛋白错误改写为两个 GH1 结构域。
            resolved = domain_family
            output_rows.append({
                "sample_id": sample_id, "contig_id": field(cgc, "contig_id"), "gene_id": gene,
                "protein_sequence_id": gene, "protein_cluster_id": clusters.get(gene, gene),
                "domain_index": index, "domain_start": domain["domain_start"], "domain_end": domain["domain_end"],
                "cazy_class": resolved[:2], "family_id": resolved,
                "protein_primary_family": family,
                "subfamily_id": domain["domain_annotation"], "EC": ec_values(overview_row, domain),
                "annotation_evidence": evidence, "selected_method": selected_method,
                "annotation_conflict": str(conflict).lower(), "CGC_id": cgc_id,
                "CGC_gene_composition": field(cgc, "annotation"), "CAZyme_substrate": domain["domain_substrate"],
                # 【问题4：无需修改】完整保留 run_dbcan 的 PUL 同源与多数投票字段，绝不在此重投票。
                "CGC_substrate_PUL": field(substrate_row, "substrate_of_the_hit_pul"),
                "CGC_substrate_vote": field(substrate_row, "substrate_predicted_by_majority_voting_of_cazymes_in_cgc"),
                "TPM": counts.get("TPM", ""), "RPM": counts.get("RPM", ""), "RPKM": counts.get("RPKM", ""),
            })
    export_rows(output_rows, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--sample-id", required=True); parser.add_argument("--overview", required=True); parser.add_argument("--abundance", required=True); parser.add_argument("--output", required=True); parser.add_argument("--hmmer"); parser.add_argument("--dbcan-sub"); parser.add_argument("--cgc-standard"); parser.add_argument("--substrate"); parser.add_argument("--cluster-map")
    args = parser.parse_args()
    standardize(args.sample_id, Path(args.overview), Path(args.abundance), Path(args.output), *(Path(value) if value else None for value in (args.hmmer, args.dbcan_sub, args.cgc_standard, args.substrate, args.cluster_map)))
