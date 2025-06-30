#!/bin/bash
#SBATCH --job-name=clustering_analysis_reservation
#SBATCH --output=clustering_analysis_reservation_%j.out
#SBATCH --error=clustering_analysis_reservation_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --reservation=mostafavilab2_Reservation
#SBATCH --partition=reservation

# Change to working directory
cd /gpfs/data/mostafavilab/rishav.e.s.dasgupta/feature_analysis

# Setup for better progress monitoring
export PYTHONUNBUFFERED=1
export PYTHONIOENCODING=utf-8

# Force immediate output flushing
exec > >(stdbuf -oL cat)
exec 2> >(stdbuf -eL cat >&2)

echo "=== CLUSTERING ANALYSIS ON RESERVATION ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "============================================"

# Load required modules
echo "Loading modules..."
module load anaconda3

# Check input files exist
echo "Checking input files..."
ORIGINAL_FILE="/gpfs/data/mostafavilab/rishav.e.s.dasgupta/all_gene_features.tsv"
STRICT_FILE="selected_features_strict.tsv"

if [ ! -f "$ORIGINAL_FILE" ]; then
    echo "ERROR: Original features file not found: $ORIGINAL_FILE"
    exit 1
fi

if [ ! -f "$STRICT_FILE" ]; then
    echo "ERROR: Strict features file not found: $STRICT_FILE"
    exit 1
fi

echo "✓ Input files verified"

# Run analysis with unbuffered output for progress monitoring
echo "Starting Python analysis at $(date)..."

# Monitor progress in background
(
    while true; do
        sleep 300  # Check every 5 minutes
        echo "[$(date)] Progress check - Job still running..."
        
        # Check if any output files are being created
        if ls *.png *.pdf *.csv 2>/dev/null | grep -q "before_after"; then
            echo "[$(date)] Output files detected - analysis nearing completion!"
            break
        fi
        
        # Check memory usage
        if command -v free >/dev/null; then
            MEM_USED=$(free -h | awk '/^Mem:/ {print $3}')
            echo "[$(date)] Memory in use: $MEM_USED"
        fi
    done
) &
MONITOR_PID=$!

# Run Python with maximum output verbosity
stdbuf -oL -eL python -u before_after_clustering_analysis.py 2>&1 | while IFS= read -r line; do
    echo "[$(date +%H:%M:%S)] $line"
done

# Stop monitoring
kill $MONITOR_PID 2>/dev/null || true

# Capture exit status
PYTHON_EXIT=$?

echo "============================================"
echo "Python analysis completed at $(date)"
echo "Exit status: $PYTHON_EXIT"

# Check outputs
echo "Checking generated files..."
OUTPUT_FILES=(
    "before_after_clustering_comparison.png"
    "before_after_clustering_comparison.pdf"
    "before_after_clustering_table.png"
    "before_after_clustering_table.pdf"
    "before_after_clustering_analysis.csv"
)

SUCCESS_COUNT=0
for file in "${OUTPUT_FILES[@]}"; do
    if [ -f "$file" ]; then
        size=$(du -sh "$file" | cut -f1)
        echo "✓ $file ($size)"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "✗ Missing: $file"
    fi
done

echo "============================================"
echo "Summary:"
echo "Generated $SUCCESS_COUNT out of ${#OUTPUT_FILES[@]} expected files"
echo "Total runtime: $SECONDS seconds"
echo "End time: $(date)"

if [ $SUCCESS_COUNT -eq ${#OUTPUT_FILES[@]} ]; then
    echo "✓ Analysis completed successfully!"
    exit 0
else
    echo "⚠ Analysis completed with missing files"
    exit 1
fi
