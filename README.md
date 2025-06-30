# Mostafavi-Gene-Heritability-Analysis-Thesis
# Gene Heritability LDSC Pipeline

A computational pipeline for identifying gene-level functional annotations that contribute to complex trait heritability using LD Score Regression (LDSC).

## Overview

This pipeline implements the methodology from "Identification of Gene Features That Are Informative of Heritability" by analyzing how gene-level functional annotations contribute to trait architecture. It moves beyond variant-level analyses to examine heritability at the gene level using comprehensive biological annotations.

## Pipeline Variants

Three distinct Nextflow pipelines for different analysis scenarios:

### 1. GO Features Pipeline (`feature_pipeline_go.nf`)
- **Purpose**: Process specific Gene Ontology features by column indices
- **Input**: Comma-separated feature column numbers
- **Use case**: Targeted analysis of functional categories (1-100 features)

### 2. Batch Processing Pipeline (`feature_pipeline_batch.nf`)
- **Purpose**: Process large feature ranges in batches
- **Input**: Start and end column ranges
- **Use case**: Large-scale feature screening (100-10,000+ features)

### 3. Benchmark Pipeline (`feature_pipeline_three_go.nf`)
- **Purpose**: Process three benchmark GO features for validation
- **Input**: Fixed benchmark set (GO:0001501, GO:0045892, GO:0048732)
- **Use case**: Method validation and comparison

## Requirements

- **Nextflow** (≥21.0), **Python 3**, **LDSC**, **bedtools** (≥2.30.0), **SLURM**
- Gene feature annotations, GENCODE annotations, 1000 Genomes reference, GWAS summary statistics

## Quick Start

```bash
# GO Features Analysis (64 features)
sbatch run_batch_production_strict.sh

# Benchmark Validation (3 features)  
sbatch run_three_go_features.sh

# Large-scale Processing
nextflow run feature_pipeline_batch.nf \
    --feature_start_column 1 --feature_end_column 1000 \
    -c nextflow_production.config

# Custom GO Features
nextflow run feature_pipeline_go.nf \
    --feature_columns "149,150,151,152" \
    -c nextflow_production.config
```

## Key Files

### Pipelines
- `feature_pipeline_go.nf` - GO features analysis
- `feature_pipeline_batch.nf` - Batch processing
- `feature_pipeline_three_go.nf` - Benchmark validation

### Configuration
- `nextflow_production.config` - High-throughput settings
- `nextflow_regular_cluster.config` - Standard cluster settings
- `params_production.yml` - Parameter template

### Scripts
- `extract_feature.py` - Single feature extraction
- `extract_feature_batch.py` - Batch feature processing
- `bedtools_ldsc_direct.sh` - SNP annotation
- `multiply_annotations.sh` - Combined annotations

### Execution & Data
- `run_batch_production_strict.sh` - GO features execution
- `run_three_go_features.sh` - Benchmark execution
- `go_features.txt` - 64 GO feature indices
- `three_go_features.txt` - 3 benchmark feature indices

## Output Structure

```
output/
├── features/                 # Extracted gene windows
├── {feature_name}/
│   ├── annotations/         # Feature annotations
│   ├── ldsc/               # LD scores
│   └── heritability/       # LDSC results
└── logs/                   # Execution logs
```

## Key Results

- **13 significantly enriched GO terms** from 66 tested (p < 0.05)
- **Context-dependent coding effects** with up to 3,000-fold enrichment
- **46 GO terms with positive enrichment** vs 20 showing depletion

## Performance

- **Memory**: 8-32GB RAM per process
- **Runtime**: 2-6 hours per feature (22 chromosomes)
- **Storage**: ~50-100MB per feature

## Example Workflows

```bash
# Custom GO analysis
echo "149,150,151" > my_features.txt
COLS=$(cat my_features.txt | tr '\n' ',' | sed 's/,$//')
nextflow run feature_pipeline_go.nf --feature_columns "$COLS"

# Batch processing
nextflow run feature_pipeline_batch.nf \
    --feature_start_column 5000 --feature_end_column 6000 \
    --batch_size 100
```

## Citation

*Dasgupta, R.E.S. (2025). Identification of Gene Features That Are Informative of Heritability. Masters thesis, Biomedical Informatics.*

## License

MIT License - see LICENSE file for details.
