#!/bin/bash
# run_batch_production_strict.sh - Updated for GO features only
#SBATCH --job-name=feature_pipeline_go
#SBATCH --output=feature_pipeline_go_%j.out
#SBATCH --error=feature_pipeline_go_%j.err
#SBATCH --time=6-20:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --reservation=mostafavilab2_Reservation
#SBATCH --partition=reservation

BASE_DIR=/gpfs/data/mostafavilab/rishav.e.s.dasgupta/feature_analysis
cd $BASE_DIR

OUTPUT_BASE=${BASE_DIR}/output_strict_features

mkdir -p ${OUTPUT_BASE}/{features,ldsc_snps,logs}

CHECKPOINT_FILE=${OUTPUT_BASE}/completed_features.txt
touch $CHECKPOINT_FILE

echo "=================================================="
echo "Feature Pipeline - GO Features Only"
echo "=================================================="
echo "Start time: $(date)"
echo "Processing GO features only (64 features)"
echo "Input file: selected_features_strict.tsv"
echo "Output directory: ${OUTPUT_BASE}"
echo "=================================================="

chmod +x ${BASE_DIR}/scripts/*.sh 2>/dev/null || true
chmod +x ${BASE_DIR}/scripts/*.py 2>/dev/null || true

# Create comma-separated list of GO columns
GO_COLUMNS=$(cat go_features.txt | tr '\n' ',' | sed 's/,$//')
GO_COUNT=$(cat go_features.txt | wc -l)

PARAMS_FILE=${OUTPUT_BASE}/params_go_features.yml
cat > $PARAMS_FILE << EOL
# Parameters for GO features only (64 features)
feature_columns: "${GO_COLUMNS}"
batch_size: ${GO_COUNT}
max_parallel_features: 25
test_mode: false
checkpoint_file: $CHECKPOINT_FILE

base_dir: "${BASE_DIR}"
output_dir: "${OUTPUT_BASE}"
scripts_dir: "${BASE_DIR}/scripts"
logs_dir: "${OUTPUT_BASE}/logs"
data_dir: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta"

feature_file: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/feature_analysis/selected_features_strict.tsv"
gencode_file: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/gencode_protein_coding_all_features_fixed.bed"

ldsc_snps_dir: "${OUTPUT_BASE}/ldsc_snps"
combined_annot_dir: "${OUTPUT_BASE}/combined_annotations"
coding_annot_dir: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Coding_Only_Annotations/"
snp_folder: "/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_EUR_Phase3_plink/"
snp_file_prefix: "1000G.EUR.QC."
ldsc_path: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/ldsc"
snp_file: "/gpfs/data/mostafavilab/shared_data/LDSC_data/hm3_no_mhc_snps.txt"
frq_path: "/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_Phase3_frq/1000G.EUR.QC."
sumstats_file: "/gpfs/data/mostafavilab/GWAS/50_irnt_standing_height_both_sexes/HEIGHT.rsid.sumstats.gz"
heritability_dir: "${OUTPUT_BASE}/heritability"

static_ref_ld_paths:
  - "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Gene_Density_Annots_Final/"
  - "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Coding_Only_Annotations/"
EOL

module load nextflow

echo "Starting pipeline with 64 GO features..."
echo "Expected processing: 64 features × 22 chromosomes = 1,408 annotation jobs"

nextflow run feature_pipeline_go.nf \
    -params-file $PARAMS_FILE \
    -c nextflow_production.config \
    -with-report ${OUTPUT_BASE}/logs/report_go_features.html \
    -with-timeline ${OUTPUT_BASE}/logs/timeline_go_features.html \
    -with-dag ${OUTPUT_BASE}/logs/dag_go_features.pdf \
    -with-trace ${OUTPUT_BASE}/logs/trace_go_features.txt \
    -ansi-log false

EXIT_CODE=$?

echo "=================================================="
echo "Pipeline completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "=================================================="

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: GO features pipeline completed"
    if [ -f $CHECKPOINT_FILE ]; then
        COMPLETED=$(wc -l < $CHECKPOINT_FILE)
        echo "Completed features: $COMPLETED/64"
    fi
    
    # Check heritability results
    echo "Heritability analysis directories:"
    find ${OUTPUT_BASE} -type d -name "heritability" | head -10
    
    # Summary statistics
    echo "Final processing summary:"
    echo "- Total features processed: 64"
    echo "- Feature type: GO terms only"
    echo "- Chromosomes: 22"
    echo "- Total annotation files: $(find ${OUTPUT_BASE} -name "*_annots_chr*.tsv" | wc -l 2>/dev/null || echo "0")"
    echo "- LDSC results: $(find ${OUTPUT_BASE} -name "*.l2.ldscore.gz" | wc -l 2>/dev/null || echo "0")"
    echo "- Heritability estimates: $(find ${OUTPUT_BASE} -name "*_h2.*" | wc -l 2>/dev/null || echo "0")"
    
    touch ${OUTPUT_BASE}/PIPELINE_SUCCESS
else
    echo "ERROR: Pipeline failed. Check logs for details."
    touch ${OUTPUT_BASE}/PIPELINE_FAILED
fi

exit $EXIT_CODE
