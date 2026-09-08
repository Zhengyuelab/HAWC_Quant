from __future__ import annotations

from pathlib import Path
import shutil
import shlex
import subprocess
import time
import os
import signal
import pandas as pd

from .config import Config
from .utils import mkdir, run_cmd, require_tools
from .annotation import gff_to_gtf, combine_files, extract_gene_protein_map, gene_lengths_from_gtf
from .matrix import htseq_to_matrix, featurecounts_to_matrix, read_matrix, write_matrix, cpm, tpm
from .orthologs import run_orthofinder, extract_species_unique_genes, normalize_search_backend
from .normalization import run_count_calibration


class Pipeline:
    def __init__(self, cfg: Config, logger):
        self.cfg = cfg
        self.log = logger
        self.outdir = mkdir(cfg.outdir)
        self.refs_dir = mkdir(self.outdir / "00-references")
        self.logs_dir = mkdir(self.outdir / "logs")

    def prepare_references(self) -> dict[str, Path]:
        refs = self.cfg.raw["references"]
        target_gff = Path(refs["target_gff"]).expanduser()
        reference_gff = Path(refs["reference_gff"]).expanduser()
        target_genome = Path(refs["target_genome"]).expanduser()
        reference_genome = Path(refs["reference_genome"]).expanduser()
        target_gtf = self.refs_dir / f"{refs['target_name']}.gtf"
        reference_gtf = self.refs_dir / f"{refs['reference_name']}.gtf"
        mixed_gtf = self.refs_dir / "genome_mixed.gtf"
        mixed_fasta = self.refs_dir / "genome_mixed.fna"
        self.log.info("Converting GFF to GTF and merging reference files")
        gff_to_gtf(target_gff, target_gtf)
        gff_to_gtf(reference_gff, reference_gtf)
        combine_files([target_gtf, reference_gtf], mixed_gtf)
        combine_files([target_genome, reference_genome], mixed_fasta)
        gene_protein = self.refs_dir / "gene_id-protein_id.tsv"
        lengths = self.refs_dir / "gene_lengths.tsv"
        extract_gene_protein_map(mixed_gtf, gene_protein)
        gene_lengths_from_gtf(mixed_gtf, lengths)
        return {"mixed_fasta": mixed_fasta, "mixed_gtf": mixed_gtf, "gene_protein": gene_protein, "gene_lengths": lengths}

    def _fastp_cmd_for_sample(self, sample: str, fq1: Path, fq2: Path, out: Path, fp: dict) -> tuple[list[str], Path, Path, Path]:
        clean_r1 = out / f"{sample}-R1.paired.clean.fq.gz"
        clean_r2 = out / f"{sample}-R2.paired.clean.fq.gz"
        json_report = out / f"{sample}.fastp.json"
        html_report = out / f"{sample}.fastp.html"
        cmd = [
            "fastp", "-i", str(fq1), "-I", str(fq2),
            "-o", str(clean_r1),
            "-O", str(clean_r2),
            "--unpaired1", str(out / f"{sample}-R1.unpaired.clean.fq.gz"),
            "--unpaired2", str(out / f"{sample}-R2.unpaired.clean.fq.gz"),
            "--thread", str(fp.get("threads", self.cfg.threads)),
            "--html", str(html_report),
            "--json", str(json_report),
        ]
        if fp.get("detect_adapter_for_pe", True): cmd.append("--detect_adapter_for_pe")
        if fp.get("cut_front", True): cmd += ["--cut_front", "--cut_front_mean_quality", str(fp.get("cut_front_mean_quality", 5))]
        if fp.get("cut_tail", True): cmd += ["--cut_tail", "--cut_tail_mean_quality", str(fp.get("cut_tail_mean_quality", 5))]
        if fp.get("cut_right", True): cmd += ["--cut_right", "--cut_right_window_size", str(fp.get("cut_right_window_size", 4)), "--cut_right_mean_quality", str(fp.get("cut_right_mean_quality", 20))]
        cmd += ["--length_required", str(fp.get("length_required", 100))]
        extra = fp.get("extra", "")
        if extra:
            cmd += shlex.split(str(extra))
        return cmd, clean_r1, clean_r2, json_report

    def write_fastp_script(self) -> Path:
        if not self.cfg.get("fastp", "enabled", default=True):
            self.log.info("Skipping fastp script because fastp.enabled is false")
            return self.logs_dir / "run_fastp_all.sh"
        require_tools(["fastp"])
        out = mkdir(self.outdir / "01-clean")
        fp = self.cfg.raw.get("fastp", {})
        samples = self.cfg.samples.reset_index(drop=True)
        skip_existing = bool(fp.get("skip_existing", True))
        timeout = fp.get("timeout_sec", None)
        if timeout in ["", 0, "0", None]:
            timeout = None
        else:
            timeout = int(timeout)

        script = self.logs_dir / "run_fastp_all.sh"
        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            "echo \"[HAWCQuant] fastp bash mode started at $(date)\"",
            f"echo \"[HAWCQuant] total samples: {len(samples)}\"",
            "",
        ]
        for idx, row in samples.iterrows():
            sample = str(row["sample"])
            fq1 = Path(str(row["fq1"])).expanduser()
            fq2 = Path(str(row["fq2"])).expanduser()
            if not fq1.is_absolute():
                fq1 = (Path.cwd() / fq1).resolve()
            if not fq2.is_absolute():
                fq2 = (Path.cwd() / fq2).resolve()
            if not fq1.exists():
                raise FileNotFoundError(f"FASTQ R1 not found for sample {sample}: {fq1}")
            if not fq2.exists():
                raise FileNotFoundError(f"FASTQ R2 not found for sample {sample}: {fq2}")
            cmd, clean_r1, clean_r2, json_report = self._fastp_cmd_for_sample(sample, fq1, fq2, out, fp)
            sample_log = self.logs_dir / f"fastp.{sample}.log"
            cmd_q = " ".join(shlex.quote(str(x)) for x in cmd)
            if timeout:
                cmd_q = f"timeout {timeout} {cmd_q}"
            lines += [f"echo \"[{idx + 1}/{len(samples)}] sample {sample} started at $(date)\""]
            if skip_existing:
                lines += [
                    f"if [[ -s {shlex.quote(str(clean_r1))} && -s {shlex.quote(str(clean_r2))} && -s {shlex.quote(str(json_report))} ]]; then",
                    f"  echo \"[{idx + 1}/{len(samples)}] sample {sample} skipped because outputs already exist\"",
                    "else",
                    f"  echo '$ {cmd_q}' > {shlex.quote(str(sample_log))}",
                    f"  {cmd_q} >> {shlex.quote(str(sample_log))} 2>&1",
                    f"  echo \"[{idx + 1}/{len(samples)}] sample {sample} finished at $(date)\"",
                    "fi",
                    "",
                ]
            else:
                lines += [
                    f"echo '$ {cmd_q}' > {shlex.quote(str(sample_log))}",
                    f"{cmd_q} >> {shlex.quote(str(sample_log))} 2>&1",
                    f"echo \"[{idx + 1}/{len(samples)}] sample {sample} finished at $(date)\"",
                    "",
                ]
        lines += ["echo \"[HAWCQuant] fastp bash mode completed at $(date)\"", ""]
        script.write_text("\n".join(lines))
        script.chmod(0o755)
        self.log.info("Wrote standalone fastp script: %s", script)
        return script

    def run_fastp(self) -> None:
        if not self.cfg.get("fastp", "enabled", default=True):
            self.log.info("Skipping fastp because fastp.enabled is false")
            return
        fp = self.cfg.raw.get("fastp", {})
        execution_mode = str(fp.get("execution_mode", "python_watchdog")).lower()
        if execution_mode in {"bash", "script", "shell"}:
            script = self.write_fastp_script()
            self.log.info("Running fastp through standalone bash script")
            run_cmd(["bash", str(script)], self.log, log_path=self.logs_dir / "fastp.bash_runner.log", stream=False)
            return

        require_tools(["fastp"])
        out = mkdir(self.outdir / "01-clean")
        samples = self.cfg.samples.reset_index(drop=True)
        skip_existing = bool(fp.get("skip_existing", True))
        timeout = fp.get("timeout_sec", None)
        if timeout in ["", 0, "0"]:
            timeout = None
        timeout = int(timeout) if timeout is not None else None

        self.log.info("Starting fastp QC for %d samples", len(samples))
        for idx, row in samples.iterrows():
            sample = str(row["sample"])
            fq1 = Path(str(row["fq1"])).expanduser()
            fq2 = Path(str(row["fq2"])).expanduser()
            if not fq1.is_absolute():
                fq1 = (Path.cwd() / fq1).resolve()
            if not fq2.is_absolute():
                fq2 = (Path.cwd() / fq2).resolve()
            if not fq1.exists():
                raise FileNotFoundError(f"FASTQ R1 not found for sample {sample}: {fq1}")
            if not fq2.exists():
                raise FileNotFoundError(f"FASTQ R2 not found for sample {sample}: {fq2}")

            cmd, clean_r1, clean_r2, json_report = self._fastp_cmd_for_sample(sample, fq1, fq2, out, fp)
            sample_log = self.logs_dir / f"fastp.{sample}.log"

            if skip_existing and clean_r1.exists() and clean_r2.exists() and json_report.exists():
                self.log.info("[%d/%d] Skipping %s because clean FASTQ and JSON already exist", idx + 1, len(samples), sample)
                continue

            self.log.info("[%d/%d] Running fastp for sample %s", idx + 1, len(samples), sample)
            if execution_mode in {"python_watchdog", "watchdog", "python", "direct"}:
                self._run_fastp_with_watchdog(cmd, sample, sample_log, [clean_r1, clean_r2, json_report], fp, timeout)
            else:
                run_cmd(cmd, self.log, log_path=sample_log, timeout=timeout, stream=False)
            self.log.info("[%d/%d] Finished fastp for sample %s", idx + 1, len(samples), sample)

    def _process_cpu_percent(self, pid: int) -> float | None:
        try:
            out = subprocess.check_output(["ps", "-p", str(pid), "-o", "%cpu="], text=True, stderr=subprocess.DEVNULL).strip()
            return float(out) if out else None
        except Exception:
            return None

    def _kill_process_group(self, proc: subprocess.Popen) -> None:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            try:
                proc.wait(timeout=10)
                return
            except subprocess.TimeoutExpired:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    def _run_fastp_once_watchdog(self, cmd: list[str], sample: str, log_path: Path, timeout: int | None, stall_timeout: int, check_interval: int, cpu_threshold: float, min_runtime: int) -> tuple[bool, str]:
        cmd_display = " ".join(shlex.quote(str(x)) for x in cmd)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        start_time = time.time()
        last_size = -1
        last_active = time.time()

        with open(log_path, "a", buffering=1) as handle:
            handle.write("\n" + "=" * 80 + "\n")
            handle.write(f"[HAWCQuant] sample {sample} started at {time.strftime('%F %T')}\n")
            handle.write(f"$ {cmd_display}\n\n")
            proc = subprocess.Popen(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True, start_new_session=True)
            while True:
                rc = proc.poll()
                now = time.time()
                try:
                    size = log_path.stat().st_size
                except FileNotFoundError:
                    size = 0
                cpu = self._process_cpu_percent(proc.pid) if rc is None else None

                if size != last_size:
                    last_size = size
                    last_active = now
                elif cpu is not None and cpu > cpu_threshold:
                    last_active = now

                if rc is not None:
                    handle.write(f"\n[HAWCQuant] sample {sample} exited with code {rc} at {time.strftime('%F %T')}\n")
                    return rc == 0, f"exit_code={rc}"

                elapsed = now - start_time
                inactive = now - last_active
                if timeout and elapsed > timeout:
                    handle.write(f"\n[HAWCQuant] timeout after {int(elapsed)} sec; killing process pid={proc.pid}\n")
                    self._kill_process_group(proc)
                    return False, f"timeout_after_{int(elapsed)}s"
                if elapsed >= min_runtime and inactive >= stall_timeout:
                    handle.write(f"\n[HAWCQuant] stall detected: no log growth and CPU <= {cpu_threshold}% for {int(inactive)} sec; killing process pid={proc.pid}\n")
                    self._kill_process_group(proc)
                    return False, f"stall_after_{int(inactive)}s"
                time.sleep(check_interval)

    def _run_fastp_with_watchdog(self, cmd: list[str], sample: str, log_path: Path, outputs: list[Path], fp: dict, timeout: int | None) -> None:
        max_retries = int(fp.get("max_retries", 3))
        check_interval = int(fp.get("stall_check_sec", 30))
        stall_timeout = int(fp.get("stall_timeout_sec", 600))
        cpu_threshold = float(fp.get("cpu_threshold", 1.0))
        min_runtime = int(fp.get("min_runtime_sec", 120))
        cleanup = bool(fp.get("cleanup_on_retry", True))

        for attempt in range(1, max_retries + 2):
            self.log.info("fastp sample %s attempt %d/%d with watchdog", sample, attempt, max_retries + 1)
            if cleanup:
                for path in outputs:
                    path.unlink(missing_ok=True)
                for suffix in ["-R1.unpaired.clean.fq.gz", "-R2.unpaired.clean.fq.gz", ".fastp.html"]:
                    (outputs[0].parent / f"{sample}{suffix}").unlink(missing_ok=True)
            ok, reason = self._run_fastp_once_watchdog(cmd, sample, log_path, timeout, stall_timeout, check_interval, cpu_threshold, min_runtime)
            outputs_ok = all(p.exists() and p.stat().st_size > 0 for p in outputs)
            if ok and outputs_ok:
                self.log.info("fastp sample %s completed successfully on attempt %d", sample, attempt)
                return
            self.log.warning("fastp sample %s failed or stalled on attempt %d/%d (%s; outputs_ok=%s)", sample, attempt, max_retries + 1, reason, outputs_ok)
            if attempt <= max_retries:
                time.sleep(int(fp.get("retry_sleep_sec", 10)))
        raise RuntimeError(f"fastp failed for sample {sample} after {max_retries + 1} attempts. See log: {log_path}")

    def build_index(self, mixed_fasta: Path) -> None:
        require_tools(["bowtie2-build"])
        expected = mixed_fasta.with_suffix(mixed_fasta.suffix + ".1.bt2")
        if expected.exists():
            self.log.info("Bowtie2 index already exists")
            return
        run_cmd(["bowtie2-build", "--threads", str(self.cfg.threads), str(mixed_fasta), str(mixed_fasta)], self.log, stream=False, log_path=self.logs_dir / "bowtie2-build.log")

    def run_alignment(self, mixed_fasta: Path) -> None:
        require_tools(["bowtie2", "samtools"])
        self.build_index(mixed_fasta)
        out = mkdir(self.outdir / "02-bam")
        clean = self.outdir / "01-clean"
        align_cfg = self.cfg.raw.get("alignment", {})
        for _, row in self.cfg.samples.iterrows():
            sample = row["sample"]
            fq1 = clean / f"{sample}-R1.paired.clean.fq.gz" if self.cfg.get("fastp", "enabled", default=True) else Path(row["fq1"])
            fq2 = clean / f"{sample}-R2.paired.clean.fq.gz" if self.cfg.get("fastp", "enabled", default=True) else Path(row["fq2"])
            sam = out / f"{sample}.sam"
            raw_bam = out / f"{sample}.raw.bam"
            bowtie_log = self.logs_dir / f"bowtie2.{sample}.log"
            cmd = ["bowtie2", "-x", str(mixed_fasta), "-1", str(fq1), "-2", str(fq2), "-S", str(sam), "--threads", str(self.cfg.threads)]
            if align_cfg.get("score_min"):
                cmd += ["--score-min", str(align_cfg["score_min"])]
            if align_cfg.get("extra_bowtie2"):
                cmd += str(align_cfg["extra_bowtie2"]).split()
            run_cmd(cmd, self.log, log_path=bowtie_log, stream=False)
            final_name = out / f"{sample}.name_sorted.bam"
            final_coord = out / f"{sample}.coord_sorted.bam"
            if align_cfg.get("deduplicate", True):
                namesort = out / f"{sample}.namesort.bam"
                fixmate = out / f"{sample}.fixmate.bam"
                coordsort = out / f"{sample}.coordsort.bam"
                dedup = out / f"{sample}.dedup.bam"
                run_cmd(f"samtools view -@ {self.cfg.threads} -bS {sam} > {raw_bam}", self.log, shell=True, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "sort", "-@", str(self.cfg.threads), "-n", "-o", str(namesort), str(raw_bam)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "fixmate", "-m", str(namesort), str(fixmate)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "sort", "-@", str(self.cfg.threads), "-o", str(coordsort), str(fixmate)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "markdup", "-r", "-@", str(self.cfg.threads), str(coordsort), str(dedup)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                shutil.copy2(dedup, final_coord)
                run_cmd(["samtools", "index", str(final_coord)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "sort", "-@", str(self.cfg.threads), "-n", "-o", str(final_name), str(dedup)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
            else:
                run_cmd(f"samtools view -@ {self.cfg.threads} -bS {sam} | samtools sort -@ {self.cfg.threads} -o {final_coord} -", self.log, shell=True, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "index", str(final_coord)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
                run_cmd(["samtools", "sort", "-@", str(self.cfg.threads), "-n", "-o", str(final_name), str(final_coord)], self.log, stream=False, log_path=self.logs_dir / f"samtools.{sample}.log")
            sam.unlink(missing_ok=True)

    def run_htseq(self, mixed_gtf: Path) -> Path:
        require_tools(["htseq-count"])
        bam_dir = self.outdir / "02-bam"
        out = mkdir(self.outdir / "03-counts")
        counting = self.cfg.raw.get("counting", {})
        ht = counting.get("htseq", self.cfg.raw.get("htseq", {}))
        feature_type = counting.get("feature_type", ht.get("feature_type", "CDS"))
        id_attribute = counting.get("id_attribute", ht.get("id_attribute", "gene_id"))
        count_files, samples = [], []
        for _, row in self.cfg.samples.iterrows():
            sample = row["sample"]
            bam = bam_dir / f"{sample}.name_sorted.bam"
            count = out / f"{sample}.htseq.tsv"
            cmd = [
                "htseq-count", "-f", "bam", "-r", ht.get("order", "name"), "-s", ht.get("stranded", "no"),
                "-a", str(ht.get("min_aqual", 10)), "-t", feature_type,
                "-i", id_attribute, "-m", ht.get("mode", "intersection-strict"),
                str(bam), str(mixed_gtf)
            ]
            run_cmd(" ".join(shlex.quote(str(x)) for x in cmd) + f" > {shlex.quote(str(count))}", self.log, shell=True, stream=False, log_path=self.logs_dir / f"htseq.{sample}.log")
            count_files.append(count)
            samples.append(sample)
        matrix = out / "raw_counts.tsv"
        htseq_to_matrix(count_files, samples, matrix)
        return matrix

    def run_featurecounts(self, mixed_gtf: Path) -> Path:
        require_tools(["featureCounts"])
        bam_dir = self.outdir / "02-bam"
        out = mkdir(self.outdir / "03-counts")
        counting = self.cfg.raw.get("counting", {})
        fc = counting.get("featureCounts", {})
        feature_type = counting.get("feature_type", "CDS")
        id_attribute = counting.get("id_attribute", "gene_id")
        samples = list(self.cfg.samples["sample"].astype(str))
        bam_files = []
        for sample in samples:
            coord = bam_dir / f"{sample}.coord_sorted.bam"
            name = bam_dir / f"{sample}.name_sorted.bam"
            bam_files.append(coord if coord.exists() else name)
        fc_out = out / "featureCounts.txt"
        cmd = [
            "featureCounts", "-T", str(fc.get("threads", self.cfg.threads)),
            "-s", str(fc.get("stranded", 0)),
            "-t", feature_type,
            "-g", id_attribute,
            "-Q", str(fc.get("min_mapq", 10)),
            "-a", str(mixed_gtf),
            "-o", str(fc_out),
        ]
        if bool(fc.get("paired", True)):
            cmd.append("-p")
        if bool(fc.get("count_read_pairs", True)):
            cmd.append("--countReadPairs")
        if bool(fc.get("require_both_ends_mapped", True)):
            cmd.append("-B")
        if bool(fc.get("check_chimeric_fragments", True)):
            cmd.append("-C")
        extra = fc.get("extra", "")
        if extra:
            cmd += shlex.split(str(extra))
        cmd += [str(x) for x in bam_files]
        run_cmd(cmd, self.log, stream=False, log_path=self.logs_dir / "featureCounts.log")
        matrix = out / "raw_counts.tsv"
        featurecounts_to_matrix(fc_out, samples, matrix)
        return matrix

    def run_counting(self, mixed_gtf: Path) -> Path:
        tool = str(self.cfg.raw.get("counting", {}).get("tool", "htseq")).lower()
        if tool in {"htseq", "htseq-count"}:
            self.log.info("Counting genes with htseq-count")
            return self.run_htseq(mixed_gtf)
        if tool in {"featurecounts", "featurecount", "subread"}:
            self.log.info("Counting genes with featureCounts")
            return self.run_featurecounts(mixed_gtf)
        raise ValueError("counting.tool must be htseq or featureCounts")

    def run_orthologs(self, gene_protein: Path) -> dict[str, Path]:
        refs = self.cfg.raw["references"]
        species = [refs["reference_name"], refs["target_name"]]
        out = mkdir(self.outdir / "04-orthologs")
        if self.cfg.get("orthofinder", "enabled", default=True):
            backend, extra_tools = normalize_search_backend(self.cfg.get("orthofinder", "search", default="diamond"))
            require_tools(["orthofinder"] + extra_tools)
            self.log.info("Running OrthoFinder with search backend: %s", backend)
            res = run_orthofinder({refs["reference_name"]: Path(refs["reference_protein_faa"]), refs["target_name"]: Path(refs["target_protein_faa"])}, out, int(self.cfg.get("orthofinder", "threads", default=self.cfg.threads)), self.cfg.get("orthofinder", "search", default="diamond"), self.log)
            og_dir = res / "Orthogroups"
        else:
            og_dir = Path(self.cfg.get("orthofinder", "orthogroups_dir"))
        return extract_species_unique_genes(og_dir / "Orthogroups.tsv", og_dir / "Orthogroups.GeneCount.tsv", species, gene_protein, out)

    def _resolve_normalization_groups(self) -> tuple[dict[str, list[str]], dict]:
        """Resolve replicate/sample groups used by the calibration model.

        `normalization.replicate_groups` can be either an explicit mapping or
        `auto`. In auto mode, groups are inferred from `samples.tsv` using
        `normalization.group_column`, usually the `condition` column. This makes
        the same workflow compatible with both replicated designs and
        single-sample-per-group designs.
        """
        norm = dict(self.cfg.raw["normalization"])
        raw_groups = norm.get("replicate_groups", "auto")
        group_column = norm.get("group_column", "condition")
        samples_df = self.cfg.samples

        if raw_groups is None or raw_groups == "auto":
            if group_column in {"sample_prefix", "sample_base", "sample_group"}:
                def derive_sample_prefix(row):
                    sample = str(row["sample"])
                    rep = str(row.get("replicate", ""))
                    parts = sample.split("-")
                    # Original HAWC design uses sample IDs such as 2-1-1,
                    # where the last field is the replicate. For no-replicate
                    # designs such as 2-1, keep the sample ID as the group ID.
                    if len(parts) >= 3 and rep and parts[-1] == rep:
                        return "-".join(parts[:-1])
                    return sample
                tmp = samples_df.copy()
                tmp["__hawc_group__"] = tmp.apply(derive_sample_prefix, axis=1)
                groups = (
                    tmp.groupby("__hawc_group__", sort=False)["sample"]
                    .apply(lambda x: [str(v) for v in x])
                    .to_dict()
                )
            else:
                if group_column not in samples_df.columns:
                    raise ValueError(f"normalization.group_column not found in samples.tsv: {group_column}")
                groups = (
                    samples_df.groupby(group_column, sort=False)["sample"]
                    .apply(lambda x: [str(v) for v in x])
                    .to_dict()
                )
        elif isinstance(raw_groups, dict):
            groups = {str(k): [str(x) for x in v] for k, v in raw_groups.items()}
        else:
            raise ValueError("normalization.replicate_groups must be a mapping or 'auto'")

        if not groups:
            raise ValueError("No normalization groups were resolved. Check samples.tsv and normalization.replicate_groups.")

        def resolve_selection(value, default):
            if value is None:
                return list(groups.keys()) if default == "all" else default
            if isinstance(value, str):
                if value in {"auto", "all"}:
                    return list(groups.keys())
                if value in {"auto_mix", "mix"}:
                    if "library" not in samples_df.columns:
                        raise ValueError("reference_replicate_groups: auto_mix requires a 'library' column in samples.tsv")
                    if group_column in {"sample_prefix", "sample_base", "sample_group"}:
                        mix_samples = samples_df.loc[samples_df["library"].astype(str).str.lower().eq("mix"), "sample"].astype(str).tolist()
                        mix_groups = [g for g, members in groups.items() if any(s in members for s in mix_samples)]
                    else:
                        mix_groups = samples_df.loc[samples_df["library"].astype(str).str.lower().eq("mix"), group_column].astype(str).drop_duplicates().tolist()
                    return [g for g in mix_groups if g in groups]
            return value

        norm["target_replicate_groups"] = resolve_selection(norm.get("target_replicate_groups", "auto"), "all")
        norm["reference_replicate_groups"] = resolve_selection(norm.get("reference_replicate_groups", "auto_mix"), "auto_mix")

        self.log.info("Resolved normalization groups: %s", groups)
        self.log.info("Target replicate groups: %s", norm["target_replicate_groups"])
        self.log.info("Reference replicate groups: %s", norm["reference_replicate_groups"])
        return groups, norm

    def run_normalization(self, raw_matrix: Path, unique_gene_files: dict[str, Path]) -> Path:
        groups, norm_cfg = self._resolve_normalization_groups()
        outputs = run_count_calibration(raw_matrix, unique_gene_files, groups, norm_cfg, mkdir(self.outdir / "05-calibration"))
        return outputs["calibrated_counts"]

    def run_expression_units(self, calibrated_counts: Path, gene_lengths: Path) -> None:
        if not self.cfg.raw.get("expression", {}).get("tpm_output", True) and not self.cfg.raw.get("expression", {}).get("cpm_output", True):
            return
        out = mkdir(self.outdir / "06-expression")
        counts = read_matrix(calibrated_counts)
        if self.cfg.raw.get("expression", {}).get("cpm_output", True):
            write_matrix(cpm(counts), out / "calibrated_CPM.tsv")
        if self.cfg.raw.get("expression", {}).get("tpm_output", True):
            lengths = pd.read_csv(gene_lengths, sep="\t").set_index("gene_id")["length"]
            write_matrix(tpm(counts, lengths), out / "calibrated_TPM.tsv")

    def run_deseq2(self, calibrated_counts: Path) -> None:
        if not self.cfg.raw.get("deseq2", {}).get("enabled", False):
            return
        tool = self.cfg.raw["deseq2"].get("trinity_run_DE_analysis", "run_DE_analysis.pl")
        require_tools([tool])
        out = mkdir(self.outdir / "07-deseq2")
        samples = self.cfg.samples[["condition", "sample"]].drop_duplicates()
        samples_file = out / "samples_described.txt"
        samples.to_csv(samples_file, sep="\t", header=False, index=False)
        shutil.copy2(calibrated_counts, out / "raw_counts.matrix")
        run_cmd([tool, "--matrix", str(out / "raw_counts.matrix"), "--method", "DESeq2", "--samples_file", str(samples_file)], self.log, cwd=out, stream=False, log_path=self.logs_dir / "deseq2.log")

    def run_all(self) -> None:
        refs = self.prepare_references()
        self.run_fastp()
        self.run_alignment(refs["mixed_fasta"])
        raw = self.run_counting(refs["mixed_gtf"])
        unique = self.run_orthologs(refs["gene_protein"])
        calibrated = self.run_normalization(raw, unique)
        self.run_expression_units(calibrated, refs["gene_lengths"])
        self.run_deseq2(calibrated)
        self.log.info("HAWCQuant pipeline completed. Results: %s", self.outdir)
