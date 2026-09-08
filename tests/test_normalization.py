from pathlib import Path
import pandas as pd

from hawcquant.normalization import proportional_fit
from hawcquant.matrix import tpm, featurecounts_to_matrix
from hawcquant.orthologs import normalize_search_backend


def test_proportional_fit_zero_intercept():
    x = pd.Series([1, 2, 3, 4], name="x")
    y = pd.Series([2, 4, 6, 8], name="y")
    fit = proportional_fit(x, y)
    assert abs(fit["slope"] - 2.0) < 1e-8
    assert abs(fit["r2"] - 1.0) < 1e-8


def test_tpm_column_sums():
    counts = pd.DataFrame({"s1": [100, 100], "s2": [50, 150]}, index=["g1", "g2"])
    lengths = pd.Series({"g1": 1000, "g2": 2000})
    out = tpm(counts, lengths)
    assert all(abs(out.sum(axis=0) - 1_000_000) < 1e-6)


def test_featurecounts_to_matrix(tmp_path: Path):
    fc = tmp_path / "featureCounts.txt"
    fc.write_text(
        "# Program:featureCounts\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\tb.bam\n"
        "gene1\tchr1\t1\t10\t+\t10\t5\t6\n"
        "gene2\tchr1\t20\t30\t-\t11\t7\t8\n"
    )
    out = tmp_path / "raw_counts.tsv"
    mat = featurecounts_to_matrix(fc, ["s1", "s2"], out)
    assert list(mat.columns) == ["s1", "s2"]
    assert int(mat.loc["gene1", "s1"]) == 5


def test_blastp_maps_to_blast():
    backend, tools = normalize_search_backend("blastp")
    assert backend == "blast"
    assert "blastp" in tools

from hawcquant.normalization import run_count_calibration


def test_no_replicate_groups_skip_within_group_normalization(tmp_path: Path):
    raw = pd.DataFrame(
        {
            "2-1": [10, 20, 5, 10],
            "2-2": [30, 40, 5, 10],
            "2-3": [20, 40, 10, 20],
            "2-4": [60, 80, 20, 40],
        },
        index=["t1", "t2", "r1", "r2"],
    )
    raw.index.name = "gene_id"
    raw_path = tmp_path / "raw.tsv"
    raw.to_csv(raw_path, sep="\t")
    target_genes = tmp_path / "target.txt"
    ref_genes = tmp_path / "ref.txt"
    target_genes.write_text("t1\nt2\n")
    ref_genes.write_text("r1\nr2\n")
    groups = {"2-1": ["2-1"], "2-2": ["2-2"], "2-3": ["2-3"], "2-4": ["2-4"]}
    cfg = {
        "target_species": "LW13",
        "reference_species": "Ecoli",
        "unique_min_count": 1,
        "z_threshold": 0,
        "fit_intercept": False,
        "min_fit_genes": 1,
        "replicate_normalization": {"enabled": "auto", "skip_single_sample_groups": True},
        "target_replicate_groups": "auto",
        "reference_replicate_groups": ["2-3", "2-4"],
        "calibration_edges": [
            {"name": "k1", "species": "target", "source_group": "2-1", "target_group": "2-3"},
            {"name": "k2", "species": "target", "source_group": "2-2", "target_group": "2-4"},
            {"name": "k3", "species": "reference", "source_group": "2-3", "target_group": "2-4"},
        ],
        "calibration_apply": [
            {"groups": ["2-1"], "multiply_by": ["k1", "k3"]},
            {"groups": ["2-2"], "multiply_by": ["k2"]},
            {"groups": ["2-3"], "multiply_by": ["k3"]},
            {"groups": ["2-4"], "multiply_by": []},
        ],
    }
    out = run_count_calibration(raw_path, {"LW13": target_genes, "Ecoli": ref_genes}, groups, cfg, tmp_path / "out")
    assert out["calibrated_counts"].exists()
    target_status = (tmp_path / "out" / "target_replicate_scaling" / "replicate_normalization_status.json").read_text()
    assert "single-sample group" in target_status
    scaling = pd.read_json(tmp_path / "out" / "target_replicate_scaling" / "replicate_scaling.json", typ="series")
    assert len(scaling) == 0

from hawcquant.config import Config
from hawcquant.pipeline import Pipeline
from hawcquant.utils import setup_logger


def test_pipeline_auto_groups_from_sample_prefix(tmp_path: Path):
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tcondition\treplicate\tlibrary\tfq1\tfq2\n"
        "2-1-1\tc1\t1\ttarget\ta\tb\n"
        "2-1-2\tc1\t2\ttarget\ta\tb\n"
        "2-3-1\tc1\t1\tmix\ta\tb\n"
    )
    cfg = Config(
        raw={
            "project": {"outdir": str(tmp_path / "out")},
            "samples": str(samples),
            "normalization": {
                "group_column": "sample_prefix",
                "replicate_groups": "auto",
                "target_replicate_groups": "auto",
                "reference_replicate_groups": "auto_mix",
            },
        },
        path=tmp_path / "config.yaml",
    )
    pipe = Pipeline(cfg, setup_logger(verbose=False))
    groups, norm = pipe._resolve_normalization_groups()
    assert groups == {"2-1": ["2-1-1", "2-1-2"], "2-3": ["2-3-1"]}
    assert norm["reference_replicate_groups"] == ["2-3"]
