#!/bin/bash
# multiply_annotations.sh
feature_annot_file=$1
coding_annot_dir=$2
feature_name=$3
chr_num=$4
output_dir=$5

coding_annot_file="${coding_annot_dir}/${chr_num}.annot.gz"

if [ ! -f "$feature_annot_file" ] || [ ! -f "$coding_annot_file" ]; then
    echo "Error: Input file not found"
    exit 1
fi

output_file="${output_dir}/${feature_name}_coding_annots_chr${chr_num}.tsv"
temp_dir=$(mktemp -d)
trap "rm -rf $temp_dir" EXIT

gunzip -c "$coding_annot_file" > "$temp_dir/coding_annot.tsv"

head -n 1 "$feature_annot_file" > "$temp_dir/feature_header.txt"
head -n 1 "$temp_dir/coding_annot.tsv" > "$temp_dir/coding_header.txt"

coding_column=$(head -n 1 "$temp_dir/coding_annot.tsv" | awk '{print $NF}')
feature_column=$(head -n 1 "$feature_annot_file" | awk '{print $NF}')
combined_column="${feature_name}_coding"

echo -e "CHR\tBP\tSNP\tCM\t${combined_column}" > "$output_file"

paste "$feature_annot_file" "$temp_dir/coding_annot.tsv" | \
    tail -n +2 | \
    awk -v fcol="$feature_column" -v ccol="$coding_column" '
    BEGIN {FS="\t"; OFS="\t"}
    {
        feature_val = $5;
        coding_val = $NF;
        product = (feature_val == "1" && coding_val == "1") ? "1" : "0";
        print $1, $2, $3, $4, product;
    }' >> "$output_file"

rm -rf "$temp_dir"
exit 0
