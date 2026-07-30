"""dbCAN 注释冲突解析与多结构域展开（本次新增修复1--3）。

Evidence: dbCAN protocol Box 6 states the family preference when methods
disagree: dbCAN HMM > dbCAN-sub/eCAMI > DIAMOND. It also documents the
overview.txt and hmmer.out/dbcan-sub.hmm.out output roles. This module applies
that stated ordering only; it does not make a new classifier.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

# dbCAN HMM can emit CAZy numeric suffixes (GH43_4); dbCAN-sub/eCAMI can
# emit identifiers such as GH43_e149. Preserve both rather than dropping eCAMI.
FAMILY = re.compile(r"\b((?:GH|GT|PL|CE|AA|CBM)\d+(?:_[A-Za-z0-9]+)?)\b")
METHODS = ("dbCAN", "dbCAN_sub", "DIAMOND")


def field(row: dict[str, str], *names: str) -> str:
    """兼容不同 run_dbcan 版本中仅大小写/空格不同的列名。"""
    normalized = {key.strip().lower().replace(" ", "_").replace("#", ""): value for key, value in row.items()}
    for name in names:
        if name in normalized:
            return normalized[name]
    return ""


def families(value: str) -> list[str]:
    """【新增修复2】提取全部家族 token，禁止按第一个结构域截断。"""
    return FAMILY.findall(value or "")


def preferred_family(overview: dict[str, str]) -> tuple[str, str, dict[str, list[str]]]:
    """Return the protocol-preferred family and all raw method calls.

    The output includes every method's tokens for audit. A disagreement is not
    discarded: ``selected_method`` records which source overrode the others.
    """
    # 【新增修复1】Box 6 的明确顺序；该循环只决定蛋白主家族，
    # 不覆盖下游 domain_rows() 保留的独立结构域家族。
    calls = {method: families(field(overview, method.lower())) for method in METHODS}
    for method in METHODS:
        if calls[method]:
            return calls[method][0].split("_")[0], method, calls
    return "", "", calls


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def domain_rows(hmmer: Path | None, dbcan_sub: Path | None) -> dict[str, list[dict[str, str]]]:
    """Collect every HMM/dbCAN-sub domain hit per protein.

    A domain hit becomes an independent result row. No first-hit reduction is
    performed. Coordinates are retained when the installed run_dbcan output
    contains them. ``overview.txt`` remains the family conflict authority.
    """
    # 以 gene_id 聚集，但每一个命中都 append；同一蛋白因此可导出多行。
    result: dict[str, list[dict[str, str]]] = {}
    for source, method, path in (("hmmer", "dbCAN", hmmer), ("dbcan_sub", "dbCAN_sub", dbcan_sub)):
        if not path or not path.exists():
            continue
        for row in read_tsv(path):
            gene = field(row, "gene_id", "gene_id_", "protein_id", "query_id", "target_name")
            annotation = field(row, "hmm_profile", "hmm_name", "family", "annotation", "cazy_family")
            if not gene or not families(annotation):
                continue
            item = {"method": method, "source_file": source, "domain_annotation": annotation,
                    "domain_start": field(row, "ali_from", "query_start", "start"),
                    "domain_end": field(row, "ali_to", "query_end", "end"),
                    "domain_ec": field(row, "ec", "ec_number"),
                    "domain_substrate": field(row, "substrate")}
            result.setdefault(gene, []).append(item)
    return result


def ec_values(overview: dict[str, str], domain: dict[str, str]) -> str:
    """Use dbCAN-sub supplied EC only; an unmapped subfamily remains blank.

    The protocol maps eCAMI subfamilies to EC/substrate through curated dbCAN-sub
    data. It does not assert that every GH/GT suffix has an EC number.
    """
    # 【新增修复3】Box 6 的 EC 来自 eCAMI/dbCAN-sub；不能把一个亚家族的
    # EC 误赋给同一融合蛋白中的其他 HMM 域，也不能凭家族后缀虚构 EC。
    if domain.get("method") != "dbCAN_sub":
        return ""
    return domain.get("domain_ec") or field(overview, "ec", "ec_")
