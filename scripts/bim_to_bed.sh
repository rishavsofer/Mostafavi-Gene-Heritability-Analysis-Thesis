#!/bin/bash
# Convert BIM files to BED format for bedtools

if [ $# -ne 2 ]; then
    echo "Usage: $0 <bim_template> <output_dir>"
    exit 1
fi

bim_template=$1  # Template like /path/to/1000G.EUR.QC.HM3_NO_HLA_SUBSET_
output_dir=$2

mkdir -p $output_dir

# Process each chromosome
for chr in {1..22}; do
    input_bim="${bim_template}${chr}.bim"
    output_bed="${output_dir}/ldsc_snps_chr${chr}.bed"
    
    echo "Converting ${input_bim} to BED format..."
    
    # BIM format: chr SNP_id genetic_dist pos ref alt
    # BED format: chr start end name (0-based, end exclusive)
    awk -v chr=$chr '{
        # Adjust for 0-based BED format (start = pos-1, end = pos)
        print "chr"$1"\t"$4-1"\t"$4"\t"$2
    }' $input_bim > $output_bed
    
    echo "Created ${output_bed} with $(wc -l < $output_bed) SNPs"
done

echo "Conversion complete!"
