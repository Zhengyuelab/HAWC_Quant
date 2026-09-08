# HAWCQuant

HAWCQuant is a one-stop, configuration-driven workflow for quantitative transcriptome calibration in experiments that contain target-species RNA-seq libraries and target-plus-reference mixed RNA-seq libraries.

The workflow supports quality control with fastp, alignment with bowtie2, optional PCR duplicate removal with samtools, gene-level counting with either htseq-count or featureCounts, species-specific gene discovery with OrthoFinder using either DIAMOND or BLASTP-style search, slope-based count calibration, CPM/TPM-like expression normalization, and optional DESeq2 differential analysis.

## Major changes in HAWCQuant v0.2.1

- Project name changed from `eginquant` to `hawcquant`.
- Python package changed from `eginquant` to `hawcquant`.
- Main command changed from `eginquant` to `hawcquant`.
- OrthoFinder search backend now accepts `diamond`, `blastp`, `blast`, or `mmseqs`.
- Counting backend now accepts `htseq` or `featureCounts`.
- Conda environment pins `fastp=1.3.2`.
- Replicate normalization now supports both replicated and no-replicate designs.
- `replicate_groups: auto` can infer groups from `samples.tsv`; single-sample groups are skipped automatically.

## Installation

```bash
git clone https://github.com/yourname/HAWCQuant.git
cd HAWCQuant
conda env create -f environment.yml
conda activate hawcquant
pip install -e .
```

Check installation:

```bash
hawcquant --version
hawcquant --help
fastp --version
```

## Quick start

Create a template config:

```bash
hawcquant init-config -o config.yaml
```

Edit `config.yaml` and `samples.tsv`, then run:

```bash
hawcquant all -c config.yaml
```

Run stages separately:

```bash
hawcquant qc -c config.yaml
hawcquant align -c config.yaml
hawcquant count -c config.yaml
hawcquant orthologs -c config.yaml
hawcquant normalize -c config.yaml
hawcquant expression -c config.yaml
```

## Config: OrthoFinder with DIAMOND or BLASTP

For large metatranscriptomic datasets, DIAMOND is usually the faster default:

```yaml
orthofinder:
  enabled: true
  search: diamond
  threads: 16
```

To use BLASTP-style OrthoFinder search, set:

```yaml
orthofinder:
  enabled: true
  search: blastp
  threads: 16
```

Internally, `blastp` is mapped to OrthoFinder's BLAST backend. HAWCQuant checks that `blastp` and `makeblastdb` are available before running this mode.

## Config: htseq-count or featureCounts

The default is htseq-count, matching the original workflow:

```yaml
counting:
  tool: htseq
  feature_type: CDS
  id_attribute: gene_id
  htseq:
    stranded: "no"
    mode: intersection-strict
    min_aqual: 10
    order: name
```

For larger datasets, featureCounts is usually much faster:

```yaml
counting:
  tool: featureCounts
  feature_type: CDS
  id_attribute: gene_id
  featureCounts:
    stranded: 0
    paired: true
    count_read_pairs: true
    require_both_ends_mapped: true
    check_chimeric_fragments: true
    min_mapq: 10
    extra: ""
```

Strandedness for featureCounts uses numeric values:

```text
0 = unstranded
1 = stranded
2 = reversely stranded
```

## Config: replicated or no-replicate designs

HAWCQuant no longer assumes that each group has three replicates. The same normalization module supports either design.

For standard replicated designs, keep `replicate_normalization.enabled: auto`; HAWCQuant will normalize within each group that contains two or more samples:

```yaml
normalization:
  group_column: sample_prefix
  replicate_groups: auto
  replicate_normalization:
    enabled: auto
    skip_single_sample_groups: true
  target_replicate_groups: auto
  reference_replicate_groups: auto_mix
```

For designs with one sample per group, use the same settings. HAWCQuant will detect that each resolved group has only one sample and skip the within-group normalization step automatically, while still calculating calibration coefficients such as `k1`, `k2`, and `k3` from the group-level values.

A no-replicate sample table can look like this:

```text
sample  condition  replicate  library  fq1                         fq2
2-1     2-1        1          target   data/reads/2-1-R1.fq.gz      data/reads/2-1-R2.fq.gz
2-2     2-2        1          target   data/reads/2-2-R1.fq.gz      data/reads/2-2-R2.fq.gz
2-3     2-3        1          mix      data/reads/2-3-R1.fq.gz      data/reads/2-3-R2.fq.gz
2-4     2-4        1          mix      data/reads/2-4-R1.fq.gz      data/reads/2-4-R2.fq.gz
```

An example is provided at:

```bash
examples/config_no_replicates.yaml
examples/samples_no_replicates.tsv
```

To force HAWCQuant to skip all within-group replicate normalization even when repeated samples exist, set:

```yaml
normalization:
  replicate_normalization:
    enabled: false
```

## Input files

### samples.tsv

Required columns:

```text
sample, condition, replicate, library, fq1, fq2
```

Example:

```text
sample   condition   replicate   library   fq1                       fq2
2-1-1    condition1  1           target    data/reads/2-1-1-R1.fq.gz data/reads/2-1-1-R2.fq.gz
2-3-1    condition1  1           mix       data/reads/2-3-1-R1.fq.gz data/reads/2-3-1-R2.fq.gz
```

### Reference files

HAWCQuant expects:

```text
target genome FASTA
target GFF
target protein FAA
reference genome FASTA
reference GFF
reference protein FAA
```

The target and reference genomes are merged before bowtie2 indexing. The target and reference GFF files are converted into a merged GTF for read counting.

## Output structure

```text
results/
  00-references/
    genome_mixed.fna
    genome_mixed.gtf
    gene_id-protein_id.tsv
    gene_lengths.tsv
  01-clean/
    *.paired.clean.fq.gz
    *.fastp.html
    *.fastp.json
  02-bam/
    *.coord_sorted.bam
    *.name_sorted.bam
  03-counts/
    raw_counts.tsv
    *.htseq.tsv              # when counting.tool: htseq
    featureCounts.txt        # when counting.tool: featureCounts
  04-orthologs/
    *_unique_genes_gene_id.txt
  05-calibration/
    all_counts.calibrated.tsv
    calibration/calibration_coefficients.json
  06-expression/
    calibrated_CPM.tsv
    calibrated_TPM.tsv
  07-deseq2/
```

## Migration from eginquant

Old command:

```bash
eginquant all -c config.yaml
```

New command:

```bash
hawcquant all -c config.yaml
```

Old Python import:

```python
import eginquant
```

New Python import:

```python
import hawcquant
```

If you are reusing an old `config.yaml`, replace the old `htseq:` block with the new `counting:` block or run:

```bash
hawcquant init-config -o config.new.yaml
```

and merge your sample/reference paths into the new template.

## fastp watchdog mode

The default QC execution mode keeps the one-stop Python workflow but monitors each fastp process. If the process shows no log growth and CPU usage remains below the configured threshold for a sustained period, HAWCQuant kills that fastp task, removes incomplete outputs, and retries the same sample automatically.

```yaml
fastp:
  execution_mode: python_watchdog
  max_retries: 3
  stall_check_sec: 30
  stall_timeout_sec: 600
  cpu_threshold: 1.0
```

## GitHub upload

```bash
git init
git add .
git commit -m "Initial release of HAWCQuant"
git branch -M main
git remote add origin https://github.com/yourname/HAWCQuant.git
git push -u origin main
```
