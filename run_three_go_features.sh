#!/bin/bash
#SBATCH --job-name=pipeline_three_go
#SBATCH --output=pipeline_three_go_%j.out
#SBATCH --error=pipeline_three_go_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_medium

BASE_DIR=/gpfs/data/mostafavilab/rishav.e.s.dasgupta/feature_analysis
cd $BASE_DIR

OUTPUT_BASE=${BASE_DIR}/output_strict_features
LOGS_DIR=${OUTPUT_BASE}/logs_three_go
mkdir -p $LOGS_DIR

CHECKPOINT_FILE=${OUTPUT_BASE}/completed_features.txt
touch $CHECKPOINT_FILE

echo "Feature Pipeline - Three GO Features: GO:0001501, GO:0045892, GO:0048732"
echo "Start: $(date)"

chmod +x ${BASE_DIR}/scripts/*.sh 2>/dev/null || true
chmod +x ${BASE_DIR}/scripts/*.py 2>/dev/null || true

THREE_GO_COLUMNS="128,3888,4278"

PARAMS_FILE=${LOGS_DIR}/params_three_go_features.yml
cat > $PARAMS_FILE << EOL
feature_columns: "${THREE_GO_COLUMNS}"
batch_size: 3
max_parallel_features: 3
test_mode: false
checkpoint_file: $CHECKPOINT_FILE
base_dir: "${BASE_DIR}"
output_dir: "${OUTPUT_BASE}"
scripts_dir: "${BASE_DIR}/scripts"
logs_dir: "${LOGS_DIR}"
data_dir: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta"
feature_file: "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/all_gene_features.tsv"
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

nextflow run feature_pipeline_three_go.nf \
    -params-file $PARAMS_FILE \
    -c nextflow_regular_cluster.config \
    -with-report ${LOGS_DIR}/report_three_go_features.html \
    -with-timeline ${LOGS_DIR}/timeline_three_go_features.html \
    -with-dag ${LOGS_DIR}/dag_three_go_features.pdf \
    -with-trace ${LOGS_DIR}/trace_three_go_features.txt \
    -ansi-log false

EXIT_CODE=$?

echo "Pipeline completed: $(date), Exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS"
    echo "129" >> $CHECKPOINT_FILE
    echo "3889" >> $CHECKPOINT_FILE
    echo "4279" >> $CHECKPOINT_FILE
    sort -u $CHECKPOINT_FILE -o $CHECKPOINT_FILE
    
    for feature in "GO_0001502" "GO_0045893" "GO_0048733"; do
        if [ -d "${OUTPUT_BASE}/${feature}" ]; then
            annot_count=$(find "${OUTPUT_BASE}/${feature}/annotations" -name "*_annots_chr*.tsv" 2>/dev/null | wc -l)
            ldsc_count=$(find "${OUTPUT_BASE}/${feature}/ldsc" -name "*.l2.ldscore.gz" 2>/dev/null | wc -l)
            h2_count=$(find "${OUTPUT_BASE}/${feature}/heritability" -name "*_h2.*" 2>/dev/null | wc -l)
            echo "$feature: ${annot_count}/22 annotations, ${ldsc_count}/22 LDSC, ${h2_count} heritability"
        else
            echo "$feature: Directory not found"
        fi
    done
    
    touch ${LOGS_DIR}/PIPELINE_SUCCESS
else
    echo "ERROR: Pipeline failed"
    touch ${LOGS_DIR}/PIPELINE_FAILED
fi

exit $EXIT_CODE
