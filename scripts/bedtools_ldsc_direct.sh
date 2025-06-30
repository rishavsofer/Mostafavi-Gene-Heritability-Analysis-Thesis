#!/bin/bash
# bedtools_ldsc_direct.sh - Direct annotation using LDSC SNPs

gene_windows_file=$1
ldsc_snps_file=$2
feature_name=$3
chr_num=$4
output_dir=$5

# Ensure consistent chromosome naming
chr_name="chr${chr_num}"

module load bedtools/2.30.0

output_file="${output_dir}/${feature_name}_annots_chr${chr_num}.tsv"

echo "Running bedtools intersect with LDSC SNPs for chromosome ${chr_num}"
echo "Gene windows: ${gene_windows_file}"
echo "LDSC SNPs file: ${ldsc_snps_file}"
echo "Using chromosome identifier: ${chr_name}"

if [ ! -f "$gene_windows_file" ]; then
    echo "Error: Gene windows file not found"
    exit 1
fi

if [ ! -f "$ldsc_snps_file" ]; then
    echo "Error: LDSC SNPs file not found"
    exit 1
fi

# Create annotation file header
echo -e "CHR\tBP\tSNP\tCM\t${feature_name}" > $output_file

# Run bedtools intersection
bedtools intersect -a ${ldsc_snps_file} -b ${gene_windows_file} -c > ${output_dir}/tmp_intersect_${chr_num}.bed

# Format the output to match LDSC requirements
# Input: chr1 752565 752566 rs3094315 1
# Output: 1 752566 rs3094315 0 1
awk -v chr=$chr_num 'BEGIN{FS="\t"; OFS="\t"}{
    # Strip "chr" prefix 
    chrnum = substr($1, 4)
    # Format as: CHR BP SNP CM Annotation
    print chrnum, $3, $4, "0", ($5 > 0 ? "1" : "0")
}' ${output_dir}/tmp_intersect_${chr_num}.bed >> $output_file

rm ${output_dir}/tmp_intersect_${chr_num}.bed

echo "SNP annotation file created: ${output_file}"
echo "Total SNPs: $(wc -l < $output_file)"
echo "SNPs with annotation=1: $(grep -c $'\t1$' $output_file)"
