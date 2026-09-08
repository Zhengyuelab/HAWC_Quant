# Migration from eginquant to HAWCQuant

## Command change

```bash
# old
eginquant all -c config.yaml

# new
hawcquant all -c config.yaml
```

## Package change

```python
# old
import eginquant

# new
import hawcquant
```

## Configuration changes

The old top-level `htseq:` block has been replaced by `counting:`:

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

For featureCounts:

```yaml
counting:
  tool: featureCounts
  feature_type: CDS
  id_attribute: gene_id
  featureCounts:
    stranded: 0
    paired: true
    count_read_pairs: true
```

OrthoFinder now supports:

```yaml
orthofinder:
  search: diamond   # or blastp, blast, mmseqs
```
