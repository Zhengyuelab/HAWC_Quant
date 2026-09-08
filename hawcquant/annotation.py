from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable
import pandas as pd


def parse_gff_attributes(attr: str) -> dict[str, str]:
    out = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        elif " " in item:
            k, v = item.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        out[k] = v
    return out


def gff_to_gtf(gff: Path, out: Path, feature_type: str = "CDS", id_attribute: str = "locus_tag") -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(gff) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != feature_type:
                continue
            attrs = parse_gff_attributes(parts[8])
            gene_id = attrs.get(id_attribute) or attrs.get("locus_tag") or attrs.get("gene") or attrs.get("ID")
            if not gene_id:
                continue
            protein_id = attrs.get("protein_id", gene_id)
            product = attrs.get("product", "unknown").replace('"', "'")
            parts[8] = f'gene_id "{gene_id}"; transcript_id "{gene_id}"; protein_id "{protein_id}"; product "{product}";'
            fout.write("\t".join(parts) + "\n")


def combine_files(files: Iterable[Path], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fout:
        for fp in files:
            with open(fp) as fin:
                for line in fin:
                    fout.write(line)


def extract_gene_protein_map(gtf: Path, out: Path) -> pd.DataFrame:
    rows = []
    pat_gene = re.compile(r'gene_id[ =]"?([^";]+)')
    pat_prot = re.compile(r'protein_id[ =]"?([^";]+)')
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            gm = pat_gene.search(parts[8])
            pm = pat_prot.search(parts[8])
            if gm:
                gene = gm.group(1).strip()
                protein = pm.group(1).strip() if pm else gene
                rows.append((gene, protein))
    df = pd.DataFrame(sorted(set(rows)), columns=["gene_id", "protein_id"])
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    return df


def gene_lengths_from_gtf(gtf: Path, out: Path) -> pd.DataFrame:
    rows = []
    pat_gene = re.compile(r'gene_id[ =]"?([^";]+)')
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            gm = pat_gene.search(parts[8])
            if gm:
                start, end = int(parts[3]), int(parts[4])
                rows.append((gm.group(1), max(1, end - start + 1)))
    df = pd.DataFrame(rows, columns=["gene_id", "length"])
    if df.empty:
        raise ValueError(f"No gene lengths parsed from {gtf}")
    df = df.groupby("gene_id", as_index=False)["length"].sum()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    return df
