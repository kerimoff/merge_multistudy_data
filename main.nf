#!/usr/bin/env nextflow
/*
========================================================================================
                         kerimoff/qtlmap
========================================================================================
 eQTL-Catalogue/qtlmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/eQTL-Catalogue/qtlmap
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     eQTL-Catalogue/qtlmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf\
     -profile tartu_hpc\
     --studyFile testdata/multi_test.tsv\
     --vcf_has_R2_field FALSE\
     --varid_rsid_map_file testdata/varid_rsid_map.tsv.gz\
     --n_batches 25

    Mandatory arguments:
      --studyFile                   Path to the TSV file containing pipeline inputs (VCF, expression matrix, metadata)

    Executions:
      --n_batches                   Number of parallel batches to run QTL Mapping per sample, must exceed the number of chromosomes (default: 400)
      --vcf_has_R2_field            Does the genotype VCF file contain R2 value in the INFO field? (default: true)
      --run_permutation             Calculate permuation p-values for each phenotype group (group_id in the phenotype metadata file) (default: false)
      --run_nominal                 Calculate nominal p-values for each phenotype group (group_id in the phenotype metadata file) (default: true)
      --n_permutations              Number of permutations to be performed per gene when run_permutation = true (default: 1000)

    QTL mapping:
      --cis_window                  The window where to search for associated variants around the phenotype (default: 1000000)
      --n_geno_pcs                  Number of genotype matrix principal components included as covariates in QTL analysis (default: 6).
      --n_pheno_pcs                 Number of phenotype matrix principal components included as covariates in QTL analysis (default: 6).
      --mincisvariant               Minimal numner of variants needed to be found in cis_window of each phenotype (default: 5)
      --covariates                  Comma-separated list of additional covariates included in the analysis (e.g. sex, age, batch). Columns with the exact same names should exist in the sample metadata file. 

    Fine mapping (SuSiE)
      --run_susie                   Perform eQTL fine mapping with SuSiE
      --vcf_genotype_field          Field in the VCF file that is used to construct the dosage matrix. Valid options are GT and DS (default: GT). 

    Format results:
      --reformat_sumstats          Add rsid and median TPM columns to the nominal summary statistics files and perform additional formatting to make the files compatible with the eQTL Catalogue (default: true)
      --varid_rsid_map_file         TSV file mapping variant ids in CHR_POS_REF_ALT format to rsids from dbSNP.

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

eQTL-Catalogue/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'eQTL-Catalogue/qtlmap'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Study file']           = params.studyFile
summary['Working dir']          = workflow.workDir
summary['Container Engine']     = workflow.containerEngine
summary['Current home']         = "$HOME"
summary['Current user']         = "$USER"
summary['Current path']         = "$PWD"
summary['Working dir']          = workflow.workDir
summary['Output dir']           = params.outdir
summary['Script dir']           = workflow.projectDir
summary['Config Profile']       = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']        = params.awsregion
   summary['AWS Queue']         = params.awsqueue
}
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

/*
 * Create a channel for input files
 */ 


Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> file(row.vcf)}
    .set { index_vcf_ch }
    
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> file(row.sample_metadata)}
    .unique()
    .collect()
    .set { sample_metadata_ch }

Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.study_name, file(row.feature_counts) ]}  
    .unique{ it[1] }
    .view()
    .set { feature_counts_ch }


include { extract_samples_from_vcf; index_vcf; merge_vcf; filter_vcf } from './modules/vcf_manupilations'
include { select_samples; rename_file } from './modules/metadata_manupilations'

workflow {
  print("in_workflow")
  //rename the feature_counts file so that they are unique
  rename_file(feature_counts_ch)

  // merge and filter vcfs
  index_vcf(index_vcf_ch)
  merge_vcf(index_vcf.out[0].collect(), index_vcf.out[1].collect())
  filter_vcf(merge_vcf.out)

  // prepare metadata of selected cell_type(qtl_group)
  select_samples(sample_metadata_ch.collect(), rename_file.out.collect())
  extract_samples_from_vcf(filter_vcf.out, select_samples.out.genotype_ids)
}

/*
 * Completion message
 */
workflow.onComplete {
    if(!workflow.success){
      log.info "Pipeline FAILED: $workflow.runName"
    } 
    else {
      log.info "Pipeline Successful: $workflow.runName"
    }

}
