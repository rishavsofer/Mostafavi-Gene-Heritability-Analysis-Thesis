#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.base_dir = "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/feature_analysis"
params.output_dir = "${params.base_dir}/output_strict_features"
params.scripts_dir = "${params.base_dir}/scripts"
params.logs_dir = "${params.base_dir}/logs_three_go"
params.data_dir = "/gpfs/data/mostafavilab/rishav.e.s.dasgupta"
params.ldsc_snps_dir = "${params.output_dir}/ldsc_snps"
params.combined_annot_dir = "${params.output_dir}/combined_annotations"
params.coding_annot_dir = "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Coding_Only_Annotations/"

params.feature_columns = ""  // Comma-separated list
params.batch_size = 3
params.max_parallel_features = 3
params.test_mode = false
params.checkpoint_file = "${params.output_dir}/completed_features.txt"

params.feature_file = "${params.data_dir}/all_gene_features.tsv"
params.gencode_file = "${params.data_dir}/gencode_protein_coding_all_features_fixed.bed"

params.snp_folder = "/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_EUR_Phase3_plink/"
params.snp_file_prefix = "1000G.EUR.QC."
params.snp_file_suffix = ".bim"
params.ldsc_path = "${params.data_dir}/ldsc"
params.snp_file = "/gpfs/data/mostafavilab/shared_data/LDSC_data/hm3_no_mhc_snps.txt"
params.weights_path = "${params.ldsc_path}/"
params.frq_path = "/gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_Phase3_frq/1000G.EUR.QC."
params.sumstats_file = "/gpfs/data/mostafavilab/GWAS/50_irnt_standing_height_both_sexes/HEIGHT.rsid.sumstats.gz"
params.heritability_dir = "${params.output_dir}/heritability"

params.static_ref_ld_paths = [
    "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Gene_Density_Annots_Final/",
    "/gpfs/data/mostafavilab/rishav.e.s.dasgupta/Coding_Only_Annotations/"
]

process convert_bim_to_bed {
    publishDir "${params.ldsc_snps_dir}", mode: 'copy'
    
    output:
    path "*.bed"
    
    script:
    bim_template = "${params.snp_folder}${params.snp_file_prefix}"
    """
    if [ ! -f "${params.ldsc_snps_dir}/ldsc_snps_chr1.bed" ]; then
        bash ${params.scripts_dir}/bim_to_bed.sh ${bim_template} ./
    else
        # Files already exist, just link them
        ln -s ${params.ldsc_snps_dir}/*.bed ./
    fi
    """
}

process extract_feature {
    tag "col_${feature_column}"
    publishDir "${params.output_dir}/features", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    val feature_column
    
    output:
    tuple path("*_gene_windows.bed"), env(FEATURE_NAME), val(feature_column)
    
    script:
    """
    module load anaconda3
    
    python3 ${params.scripts_dir}/extract_feature.py \
        ${params.feature_file} \
        ${feature_column} \
        ${params.gencode_file} \
        .
    
    if [ -f *_gene_windows.bed ]; then
        FEATURE_NAME=\$(basename *_feature.tsv | sed 's/_feature.tsv//')
    else
        echo "Error: Feature extraction failed for column ${feature_column}."
        exit 1
    fi
    """
}

process create_annotation_direct {
    tag "${feature_name}_chr${chr_num}"
    publishDir "${params.output_dir}/${feature_name}/annotations", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(chr_num), path(gene_windows_file), val(feature_name), path(ldsc_snps)
    
    output:
    tuple val(chr_num), path("${feature_name}_annots_chr${chr_num}.tsv"), val(feature_name)
    
    script:
    """
    bash ${params.scripts_dir}/bedtools_ldsc_direct.sh \
        ${gene_windows_file} \
        ${ldsc_snps} \
        ${feature_name} \
        ${chr_num} \
        .
    """
}

process multiply_annotations {
    tag "${feature_name}_chr${chr_num}"
    publishDir "${params.output_dir}/${feature_name}/combined_annotations", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(chr_num), path(annot_file), val(feature_name)
    
    output:
    tuple val(chr_num), path("${feature_name}_coding_annots_chr${chr_num}.tsv"), val("${feature_name}_coding")
    
    script:
    """
    bash ${params.scripts_dir}/multiply_annotations.sh \
        ${annot_file} \
        ${params.coding_annot_dir} \
        ${feature_name} \
        ${chr_num} \
        .
    """
}

process compute_ldsc {
    tag "${feature_name}_chr${chr_num}"
    publishDir "${params.output_dir}/${feature_name}/ldsc", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(chr_num), path(annot_file), val(feature_name)
    
    output:
    tuple val(chr_num), path("${chr_num}.l2.ldscore.gz"), path("${chr_num}.l2.M"), path("${chr_num}.l2.M_5_50"), path("${chr_num}.annot.gz"), val(feature_name)
    
    script:
    bfile_path = "${params.snp_folder}${params.snp_file_prefix}${chr_num}"
    """
    # Copy annotation file with proper naming
    cp ${annot_file} ${chr_num}.annot
    gzip -c ${annot_file} > ${chr_num}.annot.gz
    
    module load anaconda3
    source activate ldsc
    
    cd ${params.ldsc_path}
    
    # Use gzipped annotation file to preserve column names
    python2 ldsc.py --l2 --ld-wind-cm 1 \
        --bfile ${bfile_path} \
        --annot \$OLDPWD/${chr_num}.annot.gz \
        --out \$OLDPWD/${chr_num} \
        --print-snps ${params.snp_file}
    """
}

process compute_ldsc_combined {
    tag "${feature_name}_chr${chr_num}"
    publishDir "${params.output_dir}/${feature_name}/combined_annotations/ldsc", mode: 'copy'
    errorStrategy 'ignore'
    
    input:
    tuple val(chr_num), path(annot_file), val(feature_name)
    
    output:
    tuple val(chr_num), path("${chr_num}.l2.ldscore.gz"), path("${chr_num}.l2.M"), path("${chr_num}.l2.M_5_50"), path("${chr_num}.annot.gz"), val(feature_name)
    
    script:
    bfile_path = "${params.snp_folder}${params.snp_file_prefix}${chr_num}"
    """
    # Copy annotation file with proper naming
    cp ${annot_file} ${chr_num}.annot
    gzip -c ${annot_file} > ${chr_num}.annot.gz
    
    module load anaconda3
    source activate ldsc
    
    cd ${params.ldsc_path}
    
    # Use gzipped annotation file to preserve column names
    python2 ldsc.py --l2 --ld-wind-cm 1 \
        --bfile ${bfile_path} \
        --annot \$OLDPWD/${chr_num}.annot.gz \
        --out \$OLDPWD/${chr_num} \
        --print-snps ${params.snp_file}
    """
}

process compute_heritability {
    tag "${feature_name}"
    publishDir "${params.output_dir}/${feature_name}/heritability", mode: 'copy'
    
    input:
    tuple val(feature_name), val(chr_list)
    
    output:
    path("${feature_name}_h2.*")
    
    script:
    feature_ref_paths = [
        "${params.output_dir}/${feature_name}/ldsc/",
        "${params.output_dir}/${feature_name}_coding/combined_annotations/ldsc/"
    ]
    all_ref_paths = (feature_ref_paths + params.static_ref_ld_paths).join(',')
    
    """
    module load anaconda3
    source activate ldsc
    cd ${params.ldsc_path}
    
    python2 ldsc.py \
        --h2 ${params.sumstats_file} \
        --ref-ld-chr ${all_ref_paths} \
        --w-ld-chr /gpfs/data/mostafavilab/shared_data/LDSC_data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --overlap-annot \
        --frqfile-chr ${params.frq_path} \
        --print-cov \
        --print-coefficients \
        --out \$OLDPWD/${feature_name}_h2
    """
}

workflow {
    // Parse feature columns from parameter
    feature_columns = Channel.from(params.feature_columns.split(',').collect { it.toInteger() })
    
    if (params.test_mode) {
        feature_columns = feature_columns.take(params.batch_size)
    }
    
    // Check if LDSC SNPs already exist, if not create them
    ldsc_snps = convert_bim_to_bed()
    
    feature_info = extract_feature(feature_columns)
    
    chromosomes = Channel.from(1..22)
    
    ldsc_snps_by_chr = ldsc_snps
        .flatten()
        .map { file -> 
            def chrMatch = (file.name =~ /chr(\d+)/)
            def chr = chrMatch ? chrMatch[0][1].toInteger() : 0
            return [chr, file]
        }
    
    feature_with_chr = feature_info
        .combine(chromosomes)
        .map { gene_windows, feature_name, feature_col, chr -> 
            return [chr, gene_windows, feature_name, feature_col] 
        }
    
    feature_with_snps = feature_with_chr
        .combine(ldsc_snps_by_chr, by: 0)
        .map { chr, gene_windows, feature_name, feature_col, ldsc_snp ->
            return [chr, gene_windows, feature_name, ldsc_snp]
        }
    
    annotations = create_annotation_direct(feature_with_snps)
    
    combined_annotations = multiply_annotations(annotations)
    
    ldsc_results = compute_ldsc(annotations)
    ldsc_combined_results = compute_ldsc_combined(combined_annotations)
    
    // GROUP BOTH LDSC RESULTS BY FEATURE
    ldsc_regular_grouped = ldsc_results
        .map { chr_num, l2_file, m_file, m550_file, annot_gz, feature_name ->
            [feature_name, chr_num]
        }
        .groupTuple()
        .map { feature_name, chr_list ->
            [feature_name, "regular_complete"]
        }
    
    ldsc_combined_grouped = ldsc_combined_results
        .map { chr_num, l2_file, m_file, m550_file, annot_gz, feature_name ->
            // Remove "_coding" suffix to get original feature name
            def original_feature_name = feature_name.replaceAll("_coding\$", "")
            [original_feature_name, chr_num]
        }
        .groupTuple()
        .map { feature_name, chr_list ->
            [feature_name, "combined_complete"]
        }
    
    // WAIT FOR BOTH TO COMPLETE BEFORE HERITABILITY
    both_ldsc_complete = ldsc_regular_grouped
        .join(ldsc_combined_grouped)
        .map { feature_name, regular_status, combined_status ->
            [feature_name, "both_complete"]
        }

    heritability = compute_heritability(both_ldsc_complete)
}
