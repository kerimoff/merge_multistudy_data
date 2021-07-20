nextflow.enable.dsl=2

process select_samples{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'
    tag "select_samples"
    publishDir "${params.outdir}/metadata_and_fc_merged", mode: 'copy'

    input:
    path(metadata_files)
    path(feature_counts_files)

    output:
    path "merged_metadata.tsv"
    path "joined_feature_counts.tsv"
    path "genotype_ids.tsv", emit: genotype_ids

    script:
    """
    mkdir metadata_files_dir
    mkdir feature_counts_files_dir
    mv ${metadata_files.join(' ')} metadata_files_dir/
    mv ${feature_counts_files.join(' ')} feature_counts_files_dir/

    Rscript ${projectDir}/bin/merge_fc_and_metadata.R \
        -s ./metadata_files_dir\
        -e ./feature_counts_files_dir\
        -c ${params.cell_type}
    """
}
