from __future__ import annotations

from pathlib import Path
import json
from typing import Any
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .matrix import read_matrix, write_matrix, filter_genes_by_list_and_count, cap_high_zscores


def proportional_fit(x: pd.Series, y: pd.Series, fit_intercept: bool = False) -> dict[str, float]:
    df = pd.concat([x, y], axis=1).dropna()
    df = df[(df.iloc[:, 0] > 0) & (df.iloc[:, 1] > 0)]
    if df.empty:
        raise ValueError("No positive paired points for proportional fit")
    xv = df.iloc[:, 0].astype(float).to_numpy()
    yv = df.iloc[:, 1].astype(float).to_numpy()
    if fit_intercept:
        slope, intercept = np.polyfit(xv, yv, 1)
        pred = intercept + slope * xv
    else:
        slope = float(np.sum(xv * yv) / np.sum(xv * xv))
        intercept = 0.0
        pred = slope * xv
    sst = float(np.sum((yv - np.mean(yv)) ** 2))
    sse = float(np.sum((yv - pred) ** 2))
    r2 = 1 - sse / sst if sst > 0 else 1.0
    return {"slope": float(slope), "intercept": float(intercept), "r2": float(r2), "n": int(len(df))}


def plot_fit(x: pd.Series, y: pd.Series, fit: dict[str, float], out: Path, title: str) -> None:
    df = pd.concat([x, y], axis=1).dropna()
    df = df[(df.iloc[:, 0] > 0) & (df.iloc[:, 1] > 0)]
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(df.iloc[:, 0], df.iloc[:, 1], s=10, alpha=0.7)
    xs = np.logspace(np.log10(df.iloc[:, 0].min()), np.log10(df.iloc[:, 0].max()), 100)
    ys = fit["intercept"] + fit["slope"] * xs
    ax.plot(xs, ys)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(df.columns[0])
    ax.set_ylabel(df.columns[1])
    ax.set_title(title)
    ax.text(0.05, 0.95, f"y = {fit['slope']:.4g}x + {fit['intercept']:.4g}\nR² = {fit['r2']:.4f}\nn = {fit['n']}", transform=ax.transAxes, va="top")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def group_average(mat: pd.DataFrame, samples: list[str]) -> pd.Series:
    missing = [s for s in samples if s not in mat.columns]
    if missing:
        raise ValueError(f"Group contains samples absent from matrix: {missing}")
    return mat[samples].replace(0, np.nan).mean(axis=1).dropna()


def _selected_groups(groups: dict[str, list[str]], selected_groups: list[str] | str | None) -> list[str]:
    if selected_groups is None or selected_groups == "auto" or selected_groups == "all":
        return list(groups.keys())
    if isinstance(selected_groups, str):
        selected_groups = [selected_groups]
    missing = [g for g in selected_groups if g not in groups]
    if missing:
        raise ValueError(f"Selected group(s) not found in replicate_groups: {missing}")
    return list(selected_groups)


def _replicate_norm_settings(cfg: dict[str, Any]) -> tuple[str, bool]:
    raw = cfg.get("replicate_normalization", {"enabled": "auto", "skip_single_sample_groups": True})
    if isinstance(raw, dict):
        enabled = str(raw.get("enabled", "auto")).lower()
        skip_single = bool(raw.get("skip_single_sample_groups", True))
    elif isinstance(raw, bool):
        enabled = "true" if raw else "false"
        skip_single = True
    elif raw is None:
        enabled = "auto"
        skip_single = True
    else:
        enabled = str(raw).lower()
        skip_single = True
    return enabled, skip_single


def _active_replicate_groups(groups: dict[str, list[str]], selected_groups: list[str] | str | None, cfg: dict[str, Any]) -> tuple[list[str], dict[str, str]]:
    selected = _selected_groups(groups, selected_groups)
    enabled, skip_single = _replicate_norm_settings(cfg)
    status: dict[str, str] = {}

    if enabled in {"false", "no", "off", "0", "disabled", "skip"}:
        return [], {g: "skipped: replicate_normalization disabled" for g in selected}

    active: list[str] = []
    for group in selected:
        n = len(groups[group])
        if n <= 1:
            msg = "skipped: single-sample group has no within-group replicate to normalize"
            if enabled in {"true", "yes", "on", "1", "enabled"} and not skip_single:
                raise ValueError(f"Group {group} has only {n} sample(s), but replicate normalization was forced.")
            status[group] = msg
        else:
            active.append(group)
            status[group] = f"enabled: {n} samples"
    return active, status


def compute_replicate_scaling(feature_mat: pd.DataFrame, groups: dict[str, list[str]], selected_groups: list[str], outdir: Path, fit_intercept: bool = False, min_fit_genes: int = 20) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    scaled = feature_mat.copy().astype(float)
    report: dict[str, dict[str, dict[str, float]]] = {}
    outdir.mkdir(parents=True, exist_ok=True)
    for group in selected_groups:
        samples = groups[group]
        if len(samples) <= 1:
            # The caller normally filters these out; this guard keeps the function safe.
            report[group] = {}
            continue
        base = samples[0]
        report[group] = {}
        for sample in samples[1:]:
            fit = proportional_fit(feature_mat[base].rename(base), feature_mat[sample].rename(sample), fit_intercept=fit_intercept)
            if fit["n"] < min_fit_genes:
                raise ValueError(f"Too few genes for {group} {base} vs {sample}: {fit['n']} < {min_fit_genes}")
            scaled[sample] = (scaled[sample] - fit["intercept"]) / fit["slope"]
            report[group][sample] = fit
            plot_fit(feature_mat[base].rename(base), feature_mat[sample].rename(sample), fit, outdir / f"{group}_{base}_vs_{sample}.png", f"{group}: {base} vs {sample}")
    (outdir / "replicate_scaling.json").write_text(json.dumps(report, indent=2))
    return scaled, report


def apply_sample_scaling(mat: pd.DataFrame, groups: dict[str, list[str]], scaling_report: dict[str, dict[str, dict[str, float]]]) -> pd.DataFrame:
    out = mat.copy().astype(float)
    for group, sample_fits in scaling_report.items():
        for sample, fit in sample_fits.items():
            out[sample] = (out[sample] - fit.get("intercept", 0.0)) / fit["slope"]
    return out.clip(lower=0)


def compute_calibration_edges(species_mats: dict[str, pd.DataFrame], groups: dict[str, list[str]], edges: list[dict], outdir: Path, fit_intercept: bool = False, min_fit_genes: int = 20) -> dict[str, dict[str, float]]:
    outdir.mkdir(parents=True, exist_ok=True)
    coefs: dict[str, dict[str, float]] = {}
    for edge in edges:
        species = edge["species"]
        mat = species_mats[species]
        source = edge["source_group"]
        target = edge["target_group"]
        if source not in groups or target not in groups:
            raise ValueError(f"Calibration edge {edge['name']} references unknown group(s): {source}, {target}")
        x = group_average(mat, groups[source]).rename(source)
        y = group_average(mat, groups[target]).rename(target)
        fit = proportional_fit(x, y, fit_intercept=fit_intercept)
        if fit["n"] < min_fit_genes:
            raise ValueError(f"Too few genes for calibration {edge['name']}: {fit['n']} < {min_fit_genes}")
        coefs[edge["name"]] = fit | {"species": species, "source_group": source, "target_group": target}
        plot_fit(x, y, fit, outdir / f"{edge['name']}_{source}_vs_{target}.png", f"{edge['name']}: {source} -> {target}")
    (outdir / "calibration_coefficients.json").write_text(json.dumps(coefs, indent=2))
    return coefs


def apply_calibration(mat: pd.DataFrame, groups: dict[str, list[str]], coefs: dict[str, dict[str, float]], rules: list[dict]) -> pd.DataFrame:
    out = mat.copy().astype(float)
    for rule in rules:
        multiplier = 1.0
        for name in rule.get("multiply_by", []):
            multiplier *= float(coefs[name]["slope"])
        for group in rule["groups"]:
            if group not in groups:
                raise ValueError(f"Calibration rule references unknown group: {group}")
            for sample in groups[group]:
                if sample in out.columns:
                    out[sample] *= multiplier
    return out.clip(lower=0)


def _flatten_group_samples(groups: dict[str, list[str]], selected_groups: list[str] | str | None) -> list[str]:
    selected = _selected_groups(groups, selected_groups)
    return [s for g in selected for s in groups[g]]


def _write_replicate_status(outdir: Path, stage: str, selected: list[str] | str | None, active: list[str], status: dict[str, str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": stage,
        "selected_groups": selected,
        "active_groups": active,
        "status": status,
    }
    (outdir / "replicate_normalization_status.json").write_text(json.dumps(payload, indent=2))


def run_count_calibration(raw_matrix: Path, unique_gene_files: dict[str, Path], groups: dict[str, list[str]], cfg: dict, outdir: Path) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    raw = read_matrix(raw_matrix)
    min_count = float(cfg.get("unique_min_count", 10))
    z_threshold = float(cfg.get("z_threshold", 3))
    fit_intercept = bool(cfg.get("fit_intercept", False))
    min_fit_genes = int(cfg.get("min_fit_genes", 20))
    target = cfg["target_species"]
    ref = cfg["reference_species"]
    target_group_selection = cfg.get("target_replicate_groups", list(groups))
    reference_group_selection = cfg.get("reference_replicate_groups", [])

    unique_mats: dict[str, pd.DataFrame] = {}
    for sp, gene_file in unique_gene_files.items():
        genes = [x.strip() for x in Path(gene_file).read_text().splitlines() if x.strip()]
        samples = None
        if sp == ref and reference_group_selection:
            samples = _flatten_group_samples(groups, reference_group_selection)
        sub = filter_genes_by_list_and_count(raw, genes, samples, min_count)
        sub = cap_high_zscores(sub, z_threshold)
        unique_mats[sp] = sub
        write_matrix(sub, outdir / f"{sp}_unique_counts.filtered.tsv")

    # 1. Optional within-group normalization based on target-specific genes.
    active_target_groups, target_status = _active_replicate_groups(groups, target_group_selection, cfg)
    _write_replicate_status(outdir / "target_replicate_scaling", "target", target_group_selection, active_target_groups, target_status)
    if active_target_groups:
        target_scaled_unique, target_scaling = compute_replicate_scaling(
            unique_mats[target], groups, active_target_groups, outdir / "target_replicate_scaling", fit_intercept, min_fit_genes
        )
        target_replicate_scaled_all = apply_sample_scaling(raw, groups, target_scaling)
    else:
        target_scaled_unique = unique_mats[target].copy().astype(float)
        target_scaling = {}
        (outdir / "target_replicate_scaling").mkdir(parents=True, exist_ok=True)
        (outdir / "target_replicate_scaling" / "replicate_scaling.json").write_text(json.dumps(target_scaling, indent=2))
        target_replicate_scaled_all = raw.copy().astype(float)
    write_matrix(target_scaled_unique, outdir / f"{target}_unique_counts.replicate_scaled.tsv")
    write_matrix(target_replicate_scaled_all, outdir / "all_counts.target_replicate_scaled.tsv", integer=True)

    # 2. Optional within-group normalization based on reference-specific genes.
    ref_genes = [x.strip() for x in Path(unique_gene_files[ref]).read_text().splitlines() if x.strip()]
    reference_samples = _flatten_group_samples(groups, reference_group_selection) if reference_group_selection else None
    ref_unique_from_target_scaled = filter_genes_by_list_and_count(target_replicate_scaled_all, ref_genes, reference_samples, min_count)
    ref_unique_from_target_scaled = cap_high_zscores(ref_unique_from_target_scaled, z_threshold)

    active_ref_groups, ref_status = _active_replicate_groups(groups, reference_group_selection, cfg) if reference_group_selection else ([], {})
    _write_replicate_status(outdir / "reference_replicate_scaling", "reference", reference_group_selection, active_ref_groups, ref_status)
    if active_ref_groups:
        ref_scaled_unique, ref_scaling = compute_replicate_scaling(
            ref_unique_from_target_scaled, groups, active_ref_groups, outdir / "reference_replicate_scaling", fit_intercept, min_fit_genes
        )
        fully_replicate_scaled_all = apply_sample_scaling(target_replicate_scaled_all, groups, ref_scaling)
    else:
        ref_scaled_unique = ref_unique_from_target_scaled.copy().astype(float)
        ref_scaling = {}
        (outdir / "reference_replicate_scaling").mkdir(parents=True, exist_ok=True)
        (outdir / "reference_replicate_scaling" / "replicate_scaling.json").write_text(json.dumps(ref_scaling, indent=2))
        fully_replicate_scaled_all = target_replicate_scaled_all.copy().astype(float)
    write_matrix(ref_scaled_unique, outdir / f"{ref}_unique_counts.replicate_scaled.tsv")
    write_matrix(fully_replicate_scaled_all, outdir / "all_counts.replicate_scaled.tsv", integer=True)

    species_mats = {
        "target": target_scaled_unique,
        "reference": ref_scaled_unique,
        target: target_scaled_unique,
        ref: ref_scaled_unique,
    }
    coefs = compute_calibration_edges(species_mats, groups, cfg.get("calibration_edges", []), outdir / "calibration", fit_intercept, min_fit_genes)
    calibrated = apply_calibration(fully_replicate_scaled_all, groups, coefs, cfg.get("calibration_apply", []))
    write_matrix(calibrated, outdir / "all_counts.calibrated.tsv", integer=True)

    return {
        "target_replicate_scaled": outdir / "all_counts.target_replicate_scaled.tsv",
        "replicate_scaled": outdir / "all_counts.replicate_scaled.tsv",
        "calibrated_counts": outdir / "all_counts.calibrated.tsv",
        "calibration_coefficients": outdir / "calibration/calibration_coefficients.json",
    }
