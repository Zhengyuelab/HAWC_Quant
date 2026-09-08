from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml
import pandas as pd

from .utils import pathify


@dataclass
class Config:
    raw: dict[str, Any]
    path: Path

    @property
    def outdir(self) -> Path:
        return pathify(self.raw.get("project", {}).get("outdir", "results"))

    @property
    def threads(self) -> int:
        return int(self.raw.get("project", {}).get("threads", 8))

    @property
    def samples_file(self) -> Path:
        return pathify(self.raw["samples"])

    @property
    def samples(self) -> pd.DataFrame:
        df = pd.read_csv(self.samples_file, sep="\t")
        required = {"sample", "condition", "replicate", "library", "fq1", "fq2"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"samples.tsv lacks required columns: {sorted(missing)}")
        df["sample"] = df["sample"].astype(str)
        df["condition"] = df["condition"].astype(str)
        df["library"] = df["library"].astype(str)
        return df

    def get(self, *keys: str, default: Any = None) -> Any:
        node: Any = self.raw
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node


def load_config(path: str | Path) -> Config:
    p = pathify(path)
    with open(p) as handle:
        raw = yaml.safe_load(handle)
    if not isinstance(raw, dict):
        raise ValueError("Config YAML must be a dictionary.")
    return Config(raw=raw, path=p)


def write_default_config(path: str | Path) -> None:
    template = """# HAWCQuant configuration file
project:
  outdir: results
  threads: 16

# Tab-delimited file with columns:
# sample, condition, replicate, library, fq1, fq2
# library must be target or mix
samples: examples/samples.tsv

references:
  target_name: LW13
  reference_name: Ecoli
  target_genome: data/ref/LW13/genome.fna
  reference_genome: data/ref/Ecoli/genome.fna
  target_gff: data/ref/LW13/genomic.gff
  reference_gff: data/ref/Ecoli/genomic.gff
  target_protein_faa: data/ref/LW13/protein.faa
  reference_protein_faa: data/ref/Ecoli/protein.faa

fastp:
  enabled: true
  threads: 4
  skip_existing: true
  timeout_sec: 0
  execution_mode: python_watchdog
  max_retries: 3
  stall_check_sec: 30
  stall_timeout_sec: 600
  cpu_threshold: 1.0
  min_runtime_sec: 120
  retry_sleep_sec: 10
  cleanup_on_retry: true
  detect_adapter_for_pe: true
  cut_front: true
  cut_front_mean_quality: 5
  cut_tail: true
  cut_tail_mean_quality: 5
  cut_right: true
  cut_right_window_size: 4
  cut_right_mean_quality: 20
  length_required: 100
  extra: ""

alignment:
  mapper: bowtie2
  score_min: "L,-0.1,-0.1"
  deduplicate: true
  samtools_sort_memory: 2G
  extra_bowtie2: ""

# Count tool: htseq or featureCounts
counting:
  tool: htseq
  feature_type: CDS
  id_attribute: gene_id
  htseq:
    stranded: "no"
    mode: intersection-strict
    min_aqual: 10
    order: name
  featureCounts:
    stranded: 0
    paired: true
    count_read_pairs: true
    require_both_ends_mapped: true
    check_chimeric_fragments: true
    min_mapq: 10
    extra: ""

orthofinder:
  enabled: true
  # diamond is faster and recommended for large datasets; blastp maps to OrthoFinder -S blast
  search: diamond
  threads: 16
  species_order:
    - Ecoli
    - LW13

normalization:
  target_species: LW13
  reference_species: Ecoli
  unique_min_count: 10
  z_threshold: 3
  fit_intercept: false
  min_fit_genes: 20

  # Group design. Use auto to infer groups from samples.tsv by group_column.
  # sample_prefix converts 2-1-1/2-1-2/2-1-3 to group 2-1, and keeps single samples such as 2-1 as group 2-1.
  group_column: sample_prefix
  replicate_groups: auto

  # Within-group replicate normalization.
  # auto = run only for groups with >=2 samples; single-sample groups are skipped.
  # false = skip all within-group replicate normalization.
  replicate_normalization:
    enabled: auto
    skip_single_sample_groups: true

  # auto = all resolved groups; auto_mix = groups whose samples.tsv library column is mix.
  target_replicate_groups: auto
  reference_replicate_groups: auto_mix

  calibration_edges:
    - name: k1
      species: target
      source_group: 2-1
      target_group: 2-3
    - name: k2
      species: target
      source_group: 2-2
      target_group: 2-4
    - name: k3
      species: reference
      source_group: 2-3
      target_group: 2-4
  calibration_apply:
    - groups: [2-1]
      multiply_by: [k1, k3]
    - groups: [2-2]
      multiply_by: [k2]
    - groups: [2-3]
      multiply_by: [k3]
    - groups: [2-4]
      multiply_by: []

expression:
  length_source: gtf
  tpm_output: true
  cpm_output: true

deseq2:
  enabled: false
  contrasts: examples/contrasts.tsv
  trinity_run_DE_analysis: run_DE_analysis.pl
"""
    Path(path).write_text(template)
