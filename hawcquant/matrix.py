from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np


def htseq_to_matrix(count_files: list[Path], sample_names: list[str], out: Path) -> pd.DataFrame:
    if len(count_files) != len(sample_names):
        raise ValueError("count_files and sample_names must have the same length")
    series = []
    for fp, sample in zip(count_files, sample_names):
        df = pd.read_csv(fp, sep="\t", header=None, names=["gene_id", sample], comment=None)
        df = df[~df["gene_id"].astype(str).str.startswith("__")]
        df[sample] = pd.to_numeric(df[sample], errors="coerce").fillna(0)
        series.append(df.set_index("gene_id")[sample])
    mat = pd.concat(series, axis=1).fillna(0)
    mat.index.name = "gene_id"
    out.parent.mkdir(parents=True, exist_ok=True)
    mat.astype(int).to_csv(out, sep="\t")
    return mat


def featurecounts_to_matrix(featurecounts_file: Path, sample_names: list[str], out: Path) -> pd.DataFrame:
    """Convert featureCounts output into a raw count matrix.

    featureCounts outputs columns: Geneid, Chr, Start, End, Strand, Length, bam1, bam2, ...
    The last N columns are renamed using sample_names.
    """
    df = pd.read_csv(featurecounts_file, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        raise ValueError(f"featureCounts output lacks Geneid column: {featurecounts_file}")
    if len(df.columns) < 6 + len(sample_names):
        raise ValueError("featureCounts output has fewer count columns than expected samples")
    count_cols = list(df.columns[-len(sample_names):])
    mat = df[["Geneid"] + count_cols].copy()
    mat = mat.rename(columns={"Geneid": "gene_id", **dict(zip(count_cols, sample_names))})
    mat = mat.set_index("gene_id")
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0)
    mat.index.name = "gene_id"
    out.parent.mkdir(parents=True, exist_ok=True)
    mat.astype(int).to_csv(out, sep="\t")
    return mat


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0)


def write_matrix(df: pd.DataFrame, path: Path, integer: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    if integer:
        out = out.round().astype(int)
    out.index.name = "gene_id"
    out.to_csv(path, sep="\t")


def filter_genes_by_list_and_count(mat: pd.DataFrame, genes: list[str], samples: list[str] | None, min_count: float) -> pd.DataFrame:
    keep = [g for g in genes if g in mat.index]
    sub = mat.loc[keep].copy()
    cols = samples if samples else list(sub.columns)
    missing = [c for c in cols if c not in sub.columns]
    if missing:
        raise ValueError(f"Samples not found in matrix: {missing}")
    return sub[(sub[cols] >= min_count).all(axis=1)]


def cap_high_zscores(mat: pd.DataFrame, threshold: float = 3.0) -> pd.DataFrame:
    if threshold is None or threshold <= 0:
        return mat.copy()
    out = mat.copy().astype(float)
    mean = out.mean(axis=1)
    sd = out.std(axis=1, ddof=1).replace(0, np.nan)
    z = out.sub(mean, axis=0).div(sd, axis=0)
    cap = mean + threshold * sd
    for col in out.columns:
        mask = z[col] > threshold
        out.loc[mask, col] = cap[mask]
    return out.fillna(mat)


def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return counts.div(lib, axis=1) * 1e6


def tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    common = counts.index.intersection(lengths.index)
    sub = counts.loc[common].astype(float)
    length_kb = lengths.loc[common].astype(float) / 1000.0
    rpk = sub.div(length_kb, axis=0)
    denom = rpk.sum(axis=0).replace(0, np.nan)
    return rpk.div(denom, axis=1) * 1e6
