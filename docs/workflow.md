# HAWCQuant workflow notes

HAWCQuant formalizes a quantitative transcriptome calibration design with target-only and target-plus-reference mixed RNA-seq libraries.

## Core idea

1. Use target unique genes to correct within-group replicate count differences.
2. Use reference unique genes in mixed libraries to estimate between-group scaling factors.
3. Apply the resulting coefficients to target-species counts to obtain calibrated counts.
4. Normalize calibrated counts by gene length to obtain TPM-like expression values.

## Stages

| Stage | Command |
|---|---|
| fastp quality control | `hawcquant qc` |
| bowtie2 and samtools processing | `hawcquant align` |
| htseq-count or featureCounts matrix | `hawcquant count` |
| OrthoFinder unique gene screening | `hawcquant orthologs` |
| slope-based calibration | `hawcquant normalize` |
| CPM and TPM-like expression matrix | `hawcquant expression` |
| optional DESeq2 | `hawcquant deseq2` |

## Large metatranscriptomic datasets

Use `featureCounts` for counting and `diamond` for OrthoFinder by default:

```yaml
counting:
  tool: featureCounts

orthofinder:
  search: diamond
```

Use BLASTP-style OrthoFinder if needed:

```yaml
orthofinder:
  search: blastp
```

## Replicated and no-replicate designs

HAWCQuant supports both replicated and no-replicate group designs.

```yaml
normalization:
  group_column: sample_prefix
  replicate_groups: auto
  replicate_normalization:
    enabled: auto
    skip_single_sample_groups: true
```

With `group_column: sample_prefix`, sample IDs such as `2-1-1`, `2-1-2`, and `2-1-3` are assigned to group `2-1`. If the sample ID is already a group-level sample such as `2-1`, it remains group `2-1`. In `auto` mode, within-group replicate normalization is performed only for groups containing at least two samples. Single-sample groups are carried forward unchanged and are still used for between-group calibration.
