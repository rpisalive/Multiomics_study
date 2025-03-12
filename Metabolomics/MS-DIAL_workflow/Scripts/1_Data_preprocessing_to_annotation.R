#This workflow covers the following steps:
#1 Pre-processing ()of peak table generated from MS-DIAL
#2 Naming correction of GNPS output
#3 MS1 annotation with HMDB
#Step 1 and 2 are independent, Step 2 requires submission of GNPS export from MS-DIAL to GNPS
#The functions are modified from the original one to allow processing of single ion mode data and data without QC and blank.
#The workflow was adopted from https://github.com/respiratory-immunology-lab/metabolome-lipidome-MSDIAL/blob/main/instructions_v4.7.md

setwd('...')
source('pmp_preprocess.R')
source('pmp_metab_isNamed.R')
source('gnps_format_names.R')
source('gnps_SE_names.R')
source('add_hmdb.R')
source('compare_annotations_met.R')
source('keep_annotated_met.R')
source('save_curation_table.R')
pos_output <- '...' #Set NULL if you do not have this mode ran
neg_output <- NULL #Set NULL if you do not have this mode ran
metadata <- read.csv("metadata_metabolomics.csv", row.names = 1)

metab_pos_pmp <- pmp_preprocess(pos_output, neg_output, metadata = metadata, samples_key = 'IMSCCS', intens_cols = 33:53, 
                                info_cols = 1:32, pca_group = NULL, PCA_Title = NULL, blank_name = NULL, qc_name = "QC",
                                blankFC = 5, max_perc_mv = 0.8, missingPeaksFraction = 0.8, max_rsd = 25, 
                                mv_imp_rowmax = 0.7, mv_imp_colmax = 0.7, mv_imp_method = 'knn')

metab_pos_pmp$PCA_plot # view the PCA plot
metab_pos_pmp$glog_plot # confirm glog lambda value converged at the minima
metab_pos_pmp$filtering_dimensions # view how feature/samples numbers decreased with filtering
pmp_metab_isNamed(metab_pos_pmp$glog_results)

#The following starts with GNPS output
# Import GNPS data
gnps_pos <- read_csv('...')
gnps_neg <- NULL

# Format GNPS compound names and return a joined pos and neg data.frame
gnps <- gnps_format_names(gnps_pos_df = gnps_pos, gnps_neg_df = gnps_neg)

# Add the GNPS compound names to both the glog-transformed and imputed SE objects
rowData(metab_pos_pmp$glog_results)$compound_name_gnps <- gnps_SE_names(gnps_df = gnps, 
                                                                          metab_SE = metab_pos_pmp$glog_results)
rowData(metab_pos_pmp$imputed_results)$compound_name_gnps <- gnps_SE_names(gnps_df = gnps,
                                                                             metab_SE = metab_pos_pmp$imputed_results)
hmdb_df <- readRDS('hmdb_metabolites_detect_quant_v5_20231102.rds')

# Search annotations in HMDB and add to the SE objects
metab_pos_pmp$glog_results <- add_hmdb(metab_SE = metab_pos_pmp$glog_results,
                                         hmdb = hmdb_df, mass_tol = 0.002, cores = 4)
metab_pos_pmp$imputed_results <- add_hmdb(metab_SE = metab_pos_pmp$imputed_results,
                                            hmdb = hmdb_df, mass_tol = 0.002, cores = 4)

# Prepare data.frame with alignment IDs and all four annotations and filter for at least one annotation
msdial_gnps_hmdb_glog <- compare_annotations_met(metab_pos_pmp$glog_results)
msdial_gnps_hmdb_imputed <- compare_annotations_met(metab_pos_pmp$imputed_results)

# Keep only annotated rows and generate shortname column
metab_pos_glog <- keep_annotated_met(metab_pos_pmp$glog_results)
metab_pos_imputed <- keep_annotated_met(metab_pos_pmp$imputed_results)

# Retrieve the shortname values
metab_shortnames <- rowData(metab_pos_glog)$shortname

# Get HMDB and KEGG annotations
hmdb_annotations <- rowData(metab_pos_glog)$HMDB
kegg_annotations <- rowData(metab_pos_glog)$KEGG

# Save the curation table
#save_curation_table(metab_pos_glog, 'pos_glog_curation.csv')
#save_curation_table(metab_pos_imputed, 'pos_imputed_curation.csv')

###Manual curation to pick good quality peaks

# Load quality data
metab_pos_quality <- read_csv('pos_glog_curation.csv')

#Back to MS-DIAL for manual curation for quality peaks.
#Create a column to insert tick boxes using spreadsheet software, check rows which are quality peaks (relies on personal judgement).
#Checked boxes will be converted into TRUE and can be used to filter the experimental object in the following line.
# Filter the features using the TRUE/FALSE quality_peak column
metab_pos_glog <- metab_pos_glog[metab_pos_quality$quality_peak,]

#Save .rds
saveRDS(metab_pos_glog, file = "annotated_qp.RDS") 
