"""CAZyme domain conflict arbitration and curated EC mapping.

This module implements the project rule requested for this release:
HMMER (rank 0) > Hotpep (rank 1) > DIAMOND (rank 2).  It is deliberately
separate from dbCAN's historical version-specific method preference. dbCAN-sub
is used here only as a curated EC/substrate mapping source, not as a fourth
family-call competitor.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable

FAMILY = re.compile(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+(?:_[A-Za-z0-9]+)?)\b")
# 新增修复1：全工程唯一的来源优先级定义，输出必须包含对应 source_rank。
SOURCE_RANK = {"hmmer": 0, "hotpep": 1, "diamond": 2}


def field(row: dict[str, str], *names: str) -> str:
    """兼容各 dbCAN 工具 TSV 的空格、#、大小写列名差异。"""
    normalized = {key.strip().lower().replace(" ", "_").replace("#", ""): value for key, value in row.items()}
    for name in names:
        if name in normalized:
            return normalized[name]
    return ""


def families(value: str) -> list[str]:
    """提取全部 CAZy 家族/亚家族 token；绝不按第一个命中截断。"""
    return FAMILY.findall(value or "")


def split_family(value: str) -> tuple[str, str]:
    """将 GH43_4 拆为主家族 GH43 与亚家族 GH43_4。"""
    match = FAMILY.search(value or "")
    if not match:
        return "", ""
    subfamily = match.group(1)
    return subfamily.split("_")[0], subfamily


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _coordinates(row: dict[str, str]) -> tuple[int | None, int | None]:
    """读取 domain/query 对齐区间；没有坐标的工具结果只可作保守兜底。"""
    try:
        return int(field(row, "ali_from", "query_start", "start", "from")), int(field(row, "ali_to", "query_end", "end", "to"))
    except ValueError:
        return None, None


def source_candidates(path: Path | None, source: str) -> list[dict[str, object]]:
    """把单个工具的结果转换为逐域候选，禁止用 concat 无差别并表。"""
    if source not in SOURCE_RANK or not path or not path.exists():
        return []
    candidates: list[dict[str, object]] = []
    for row in read_tsv(path):
        gene = field(row, "gene_id", "gene_id_", "protein_id", "query_id", "target_name")
        annotation = field(row, "hmm_profile", "hmm_name", "family", "annotation", "cazy_family", "cazyme")
        start, end = _coordinates(row)
        for token in families(annotation):
            family, subfamily = split_family(token)
            candidates.append({"gene_id": gene, "cazy_family": family, "cazy_subfamily": subfamily,
                               "domain_start": start, "domain_end": end, "annotation_source": source,
                               "source_rank": SOURCE_RANK[source], "raw_annotation": annotation})
    return [item for item in candidates if item["gene_id"] and item["cazy_family"]]


def interval_overlap(left: dict[str, object], right: dict[str, object]) -> float:
    """返回相交长度除以较短结构域长度，用于 >=80% 的同域判定。"""
    a_start, a_end = left["domain_start"], left["domain_end"]
    b_start, b_end = right["domain_start"], right["domain_end"]
    if None in (a_start, a_end, b_start, b_end):
        # 没有坐标时不应错误合并两个潜在独立域；仅同来源同家族才视为重复。
        return 1.0 if left["cazy_subfamily"] == right["cazy_subfamily"] else 0.0
    overlap = max(0, min(int(a_end), int(b_end)) - max(int(a_start), int(b_start)) + 1)
    shortest = min(int(a_end) - int(a_start) + 1, int(b_end) - int(b_start) + 1)
    return overlap / shortest if shortest > 0 else 0.0


def resolve_conflicts(candidates: Iterable[dict[str, object]], overlap_threshold: float = .80) -> list[dict[str, object]]:
    """按基因和域区间进行保守裁决。

    新增修复1/2/4：先按 source_rank 升序检查候选；若与已保留域的
    重叠比例 >= 阈值，则低优先级候选被丢弃。同等级或不同区间不重叠的
    最高优先级结构域均保留，因此多结构域蛋白不会被压缩为单一注释。
    """
    accepted: dict[str, list[dict[str, object]]] = {}
    order = sorted(candidates, key=lambda item: (str(item["gene_id"]), int(item["source_rank"]), int(item["domain_start"] or -1), int(item["domain_end"] or -1), str(item["cazy_subfamily"])))
    for candidate in order:
        retained = accepted.setdefault(str(candidate["gene_id"]), [])
        duplicate = next((item for item in retained if interval_overlap(candidate, item) >= overlap_threshold), None)
        if duplicate is None:
            retained.append(dict(candidate))
        # 若 duplicate 已存在，其 source_rank 不高于当前候选；当前低优先级结果不输出。
    return [item for gene in sorted(accepted) for item in accepted[gene]]


def load_ec_mapping(path: Path | None) -> dict[str, str]:
    """加载策展的 subfamily/family -> EC 表；无文件或无命中均返回空。"""
    if not path or not path.exists():
        return {}
    mapping: dict[str, str] = {}
    for row in read_tsv(path):
        key = field(row, "subfamily", "family", "cazy_subfamily", "cazy_family", "hmm_profile")
        ec = field(row, "ec", "ec_number", "ec_")
        if key and ec:
            mapping[key] = ec
    return mapping


def map_ec(candidate: dict[str, object], mapping: dict[str, str]) -> str:
    """新增修复6：先查亚家族，再降级主家族；没有策展映射则置空。"""
    return mapping.get(str(candidate["cazy_subfamily"]), mapping.get(str(candidate["cazy_family"]), ""))
