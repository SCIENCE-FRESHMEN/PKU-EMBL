"""Publication-format views corresponding to dbcan_plot [dbCAN PDF, P14--P16]."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def save(fig, output: str, dpi: int = 300) -> None:
    """Write the protocol's PDF handoff format plus a high-resolution PNG."""
    path = Path(output); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def heatmap(table: str, output: str, top: int = 20) -> None:
    """P14: top substrate TPM heatmap from fam_substrate_abund.out."""
    data = pd.read_csv(table, sep="\t")
    value_columns = data.select_dtypes("number").columns
    if not len(value_columns): raise ValueError("Expected numeric TPM abundance columns.")
    labels = data.iloc[:, 0].astype(str)
    matrix = data.loc[:, value_columns].copy(); matrix.index = labels
    matrix = matrix.loc[matrix.sum(axis=1).nlargest(top).index]
    fig, ax = plt.subplots(figsize=(10, max(4, top * .35)))
    sns.heatmap(matrix, cmap="viridis", ax=ax); ax.set_title("CAZyme substrate abundance (TPM)")
    save(fig, output)


def barplot(table: str, output: str, label: str, abundance: str, top: int = 20) -> None:
    """P15: top 20 family/subfamily/EC abundance bar plot."""
    data = pd.read_csv(table, sep="\t").sort_values(abundance, ascending=False).head(top)
    fig, ax = plt.subplots(figsize=(10, 6)); sns.barplot(data=data, y=label, x=abundance, ax=ax, color="#3274a1")
    ax.set_title("CAZyme abundance"); save(fig, output)


def cgc_structure(cgc_table: str, output: str, cgc_id: str) -> None:
    """P16 structure panel from documented cgc_standard.out columns."""
    columns = ["CGC_id", "type", "contig_id", "gene_id", "start", "end", "strand", "annotation"]
    data = pd.read_csv(cgc_table, sep="\t", names=columns, header=0)
    genes = data[data["CGC_id"].astype(str).eq(cgc_id)].sort_values("start")
    if genes.empty: raise ValueError(f"CGC not found: {cgc_id}")
    fig, ax = plt.subplots(figsize=(12, 2.5)); colors = {"CAZyme": "#55a868", "TC": "#4c72b0", "TF": "#c44e52", "STP": "#8172b3"}
    for _, gene in genes.iterrows():
        start, end = int(gene.start), int(gene.end); direction = end - start if gene.strand == "+" else start - end
        ax.arrow(start, .5, direction, 0, width=.06, head_width=.16, length_includes_head=True, color=colors.get(gene.type, "#999999"))
        ax.text((start + end) / 2, .12, str(gene.annotation), ha="center", fontsize=7, rotation=35)
    ax.set(ylim=(0, 1), yticks=[], xlabel="Contig position (bp)", title=f"CGC structure: {cgc_id}"); save(fig, output)


def pul_synteny(cgc_table: str, pul_table: str, output: str) -> None:
    """P16 schematic; full alignment requires PUL.out and dbCAN-PUL GFF as in the PDF."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 3.5), sharex=True)
    for ax, table, title in zip(axes, [cgc_table, pul_table], ["Query CGC", "Best dbCAN-PUL hit"]):
        data = pd.read_csv(table, sep="\t")
        for i, (_, row) in enumerate(data.iterrows()):
            ax.arrow(i, .5, .8, 0, width=.06, head_width=.16, length_includes_head=True)
            ax.text(i + .4, .14, str(row.iloc[-1]), ha="center", fontsize=7, rotation=35)
        ax.set(ylim=(0, 1), yticks=[], title=title)
    save(fig, output)


def main() -> None:
    parser = argparse.ArgumentParser(); sub = parser.add_subparsers(dest="mode", required=True)
    for name in ("heatmap", "bar", "cgc", "synteny"):
        p = sub.add_parser(name); p.add_argument("--output", required=True); p.add_argument("--table", required=name != "synteny")
        if name == "bar": p.add_argument("--label", required=True); p.add_argument("--abundance", required=True)
        if name == "cgc": p.add_argument("--cgc-id", required=True)
        if name == "synteny": p.add_argument("--cgc-table", required=True); p.add_argument("--pul-table", required=True)
    args = parser.parse_args()
    if args.mode == "heatmap": heatmap(args.table, args.output)
    elif args.mode == "bar": barplot(args.table, args.output, args.label, args.abundance)
    elif args.mode == "cgc": cgc_structure(args.table, args.output, args.cgc_id)
    else: pul_synteny(args.cgc_table, args.pul_table, args.output)


if __name__ == "__main__": main()
