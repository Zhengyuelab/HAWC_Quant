from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config, write_default_config
from .utils import setup_logger
from .pipeline import Pipeline
from . import __version__


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hawcquant", description="One-stop HAWC quantitative transcriptome calibration workflow")
    p.add_argument("--version", action="version", version=f"hawcquant {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init-config", help="write a template config.yaml")
    p_init.add_argument("-o", "--output", default="config.yaml")

    def add_common(sp):
        sp.add_argument("-c", "--config", required=True, help="config YAML")
        sp.add_argument("--verbose", action="store_true")

    for name in ["all", "prepare", "qc", "qc-script", "align", "count", "orthologs", "normalize", "expression", "deseq2"]:
        sp = sub.add_parser(name, help=f"run {name} stage")
        add_common(sp)
        if name in {"normalize", "expression"}:
            sp.add_argument("--raw-matrix", help="raw count matrix for standalone mode")
        if name == "normalize":
            sp.add_argument("--target-unique", help="target unique gene list")
            sp.add_argument("--reference-unique", help="reference unique gene list")
        if name == "expression":
            sp.add_argument("--calibrated-counts", help="calibrated count matrix")
    return p


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "init-config":
        write_default_config(args.output)
        print(f"Wrote template config: {args.output}")
        return

    cfg = load_config(args.config)
    logger = setup_logger(cfg.outdir / "logs" / f"hawcquant_{args.command}.log", getattr(args, "verbose", False))
    pipe = Pipeline(cfg, logger)

    if args.command == "all":
        pipe.run_all()
    elif args.command == "prepare":
        pipe.prepare_references()
    elif args.command == "qc":
        pipe.run_fastp()
    elif args.command == "qc-script":
        script = pipe.write_fastp_script()
        print(f"Wrote fastp script: {script}")
    elif args.command == "align":
        refs = pipe.prepare_references()
        pipe.run_alignment(refs["mixed_fasta"])
    elif args.command == "count":
        refs = pipe.prepare_references()
        pipe.run_counting(refs["mixed_gtf"])
    elif args.command == "orthologs":
        refs = pipe.prepare_references()
        pipe.run_orthologs(refs["gene_protein"])
    elif args.command == "normalize":
        refs = cfg.raw["references"]
        raw_matrix = Path(args.raw_matrix) if args.raw_matrix else cfg.outdir / "03-counts" / "raw_counts.tsv"
        target_unique = Path(args.target_unique) if args.target_unique else cfg.outdir / "04-orthologs" / f"{refs['target_name']}_unique_genes_gene_id.txt"
        reference_unique = Path(args.reference_unique) if args.reference_unique else cfg.outdir / "04-orthologs" / f"{refs['reference_name']}_unique_genes_gene_id.txt"
        pipe.run_normalization(raw_matrix, {refs["target_name"]: target_unique, refs["reference_name"]: reference_unique})
    elif args.command == "expression":
        refs = pipe.prepare_references()
        counts = Path(args.calibrated_counts) if args.calibrated_counts else cfg.outdir / "05-calibration" / "all_counts.calibrated.tsv"
        pipe.run_expression_units(counts, refs["gene_lengths"])
    elif args.command == "deseq2":
        counts = cfg.outdir / "05-calibration" / "all_counts.calibrated.tsv"
        pipe.run_deseq2(counts)


if __name__ == "__main__":
    main()
