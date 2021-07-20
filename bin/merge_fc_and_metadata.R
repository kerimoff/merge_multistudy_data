suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("-s", "--samplemeta_dir"), type="character", default="./metadata_files_dir",
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("-e", "--expression_matrix_dir"), type="character", default="./feature_counts_files_dir",
                        help="Expression matrix file path with gene phenotype-id in rownames and sample-is in columnnames", metavar = "type"),
  optparse::make_option(c("-o", "--outdir"), type="character", default=".",
                        help="Path to the output directory.", metavar = "type"),
  optparse::make_option(c("-c", "--cell_type"), type="character", default=NULL,
                        help="Single string value of common cell type. [default \"%default\"]", metavar = "type")
  # optparse::make_option(c("-c", "--condition"), type="character", default="naive",
  #                       help="Single string value of common condition. [default \"%default\"]", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

sample_meta_path = opt$s
expression_matrix_path = opt$e
output_dir = opt$o
common_cell_type = opt$c
# common_condition = opt$c

# assertthat::assert_that(!is.null(common_qtl_group), "QTL group should be provided.")

metadata_files <- list.files(sample_meta_path, full.names = T)
feature_counts_files <- list.files(expression_matrix_path, full.names = T)

mandatory_metadata_cols <- c("sample_id", "genotype_id", "sex", "cell_type", "condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed", "study")
merged_metadata <- base::data.frame()
for (metadata_file in metadata_files) {
  message("Processing metadata: ", metadata_file)
  metadata_tsv <- readr::read_tsv(metadata_file) %>% dplyr::filter(genotype_qc_passed, rna_qc_passed)
  assertthat::assert_that(assertthat::has_name(metadata_tsv, mandatory_metadata_cols))
  merged_metadata <- merged_metadata %>% rbind(metadata_tsv %>% select(mandatory_metadata_cols))
}

merged_metadata <- merged_metadata %>% 
  dplyr::filter(cell_type == common_cell_type) %>%
  dplyr::filter(condition == "naive") %>% 
  dplyr::mutate(qtl_group = paste0(common_cell_type, "_naive"))

message(nrow(merged_metadata))
message(merged_metadata %>% dplyr::pull(genotype_id) %>% unique() %>% length())

assertthat::assert_that(nrow(merged_metadata) == merged_metadata %>% dplyr::pull(genotype_id) %>% unique() %>% length())

base::message("writing merged metadata to: ", paste0(output_dir, "/merged_metadata.tsv"))
readr::write_tsv(merged_metadata, paste0(output_dir, "/merged_metadata.tsv"))

base::message("writing unique genotype_ids to: ", paste0(output_dir, "/genotype_ids.tsv"))
readr::write_tsv(merged_metadata %>% select(genotype_id), paste0(output_dir, "/genotype_ids.tsv"), col_names = F)


merged_feature_counts <- base::data.frame()
first_dataset <- TRUE
gene_count <- 0
for (fc_file in feature_counts_files) {
  message("Processing expression matrix: ", fc_file)
  fc_tsv <- readr::read_tsv(fc_file)
  assertthat::assert_that(assertthat::has_name(fc_tsv, "phenotype_id"))
  if (first_dataset) {
    gene_count <- nrow(fc_tsv)
    first_dataset <- FALSE
    merged_feature_counts <- fc_tsv
  } else {
    common_cols = dplyr::intersect(colnames(merged_feature_counts), colnames(fc_tsv))
    assertthat::assert_that(common_cols %>% length() == 1, msg = "There are same sample_ids in different datasets.")
    assertthat::assert_that(nrow(fc_tsv)==gene_count, msg = "gene counts are different in the datasets")
    merged_feature_counts <- merged_feature_counts %>% left_join(fc_tsv)
  }
}

base::message("writing joined expression matrix to: ", paste0(output_dir, "/joined_feature_counts.tsv"))
readr::write_tsv(merged_feature_counts, paste0(output_dir, "/joined_feature_counts.tsv"))
