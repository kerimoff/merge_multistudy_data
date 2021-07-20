nextflow.enable.dsl=2

process index_vcf{
    tag "index_vcf"
    cpus 1
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    path vcf 

    output:
    path vcf 
    path "${vcf}.csi"

    script:
    """
    bcftools index ${vcf}
    """
}

process merge_vcf {
    tag "merge_vcf"
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    file(input_vcfs) 
    file(input_vcf_indexes) 

    output:
    path "merged.vcf.gz"

    """
    bcftools merge --threads ${task.cpus} -Oz -o merged.vcf.gz ${input_vcfs.join(' ')} 
    """
}

process filter_vcf {
    tag "filter_vcf"
    container = 'quay.io/eqtlcatalogue/genimpute:v20.06.1'
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    file(vcf) 

    output:
    file ("${vcf.simpleName}_filtered.vcf.gz")

    """
    bcftools filter -i 'MAF[0] > 0.01' -Oz -o ${vcf.simpleName}_filtered.vcf.gz $vcf
    """
}

process extract_samples_from_vcf {
    tag "extract_samples"
    publishDir "${params.outdir}/vcf", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    path genotype_vcf
    path sample_names

    output:
    path ("${genotype_vcf.simpleName}_extracted.vcf.gz"), emit: vcf 

    script:
    """
    bcftools view -S $sample_names -Oz -o ${genotype_vcf.simpleName}_extracted.vcf.gz --force-samples $genotype_vcf
    """
}