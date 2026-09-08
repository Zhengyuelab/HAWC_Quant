from __future__ import annotations

from pathlib import Path
import shutil
import pandas as pd

from .utils import run_cmd


def normalize_search_backend(search: str) -> tuple[str, list[str]]:
    """Return OrthoFinder -S value and extra tools that should exist.

    Config accepts diamond, blast, blastp, mmseqs. blastp maps to OrthoFinder's blast backend.
    """
    s = str(search).strip().lower()
    if s in {"blastp", "blast"}:
        return "blast", ["blastp", "makeblastdb"]
    if s == "diamond":
        return "diamond", ["diamond"]
    if s in {"mmseqs", "mmseqs2"}:
        return "mmseqs", ["mmseqs"]
    raise ValueError("orthofinder.search must be one of: diamond, blastp, blast, mmseqs")


def run_orthofinder(faa_files: dict[str, Path], outdir: Path, threads: int = 4, search: str = "diamond", logger=None) -> Path:
    faa_dir = outdir / "faa"
    faa_dir.mkdir(parents=True, exist_ok=True)
    for name, fp in faa_files.items():
        dest = faa_dir / f"{name}.faa"
        if not dest.exists():
            shutil.copy2(fp, dest)
    backend, _ = normalize_search_backend(search)
    run_cmd(["orthofinder", "-f", str(faa_dir), "-t", str(threads), "-a", str(max(1, threads // 2)), "-S", backend], logger=logger, stream=False, log_path=outdir / "orthofinder.log")
    results = sorted(faa_dir.glob("OrthoFinder/Results_*"))
    if not results:
        results = sorted(faa_dir.glob("Results_*"))
    if not results:
        raise RuntimeError("Cannot locate OrthoFinder Results_* directory")
    return results[-1]


def _split_genes(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).replace(", ", ",").split(",") if x.strip()]


def extract_species_unique_genes(orthogroups_tsv: Path, gene_count_tsv: Path, species: list[str], gene_protein_map: Path, outdir: Path) -> dict[str, Path]:
    og = pd.read_csv(orthogroups_tsv, sep="\t")
    gc = pd.read_csv(gene_count_tsv, sep="\t")
    mapper = pd.read_csv(gene_protein_map, sep="\t")
    protein_to_gene = dict(zip(mapper["protein_id"].astype(str), mapper["gene_id"].astype(str)))
    outputs = {}
    for sp in species:
        other = [x for x in species if x != sp]
        if sp not in gc.columns:
            raise ValueError(f"Species column {sp} not found in {gene_count_tsv}. Columns: {list(gc.columns)}")
        mask = gc[sp].fillna(0).astype(int) > 0
        for o in other:
            if o not in gc.columns:
                raise ValueError(f"Species column {o} not found in {gene_count_tsv}")
            mask &= gc[o].fillna(0).astype(int) == 0
        unique_ogs = set(gc.loc[mask, "Orthogroup"].astype(str))
        proteins: list[str] = []
        for _, row in og[og["Orthogroup"].astype(str).isin(unique_ogs)].iterrows():
            proteins.extend(_split_genes(row.get(sp, "")))
        genes = sorted({protein_to_gene.get(p, p) for p in proteins})
        out = outdir / f"{sp}_unique_genes_gene_id.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(genes) + "\n")
        outputs[sp] = out
    return outputs
