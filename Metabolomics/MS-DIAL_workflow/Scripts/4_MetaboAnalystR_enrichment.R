library(MetaboAnalystR)
library(dplyr)

setwd('path_to/working_directory')
metab_SE <- readRDS('path_to/annotated_qp.RDS')
# Get HMDB and KEGG annotations
shortnames <- rowData(metab_SE)$shortname
hmdb_accession <- rowData(metab_SE)$HMDB_accession
kegg_annotations <- rowData(metab_SE)$KEGG
gnps <- rowData(metab_SE)$compound_name_gnps
chebi_ID <- rowData(metab_SE)$chebi_ID
metlin_ID <- rowData(metab_SE)$metlin_ID
pubchem_ID <- rowData(metab_SE)$pubchem_ID
annotations <- data.frame(Short_names = shortnames, HMDB_ID = hmdb_accession, KEGG_ID = kegg_annotations, chebi_ID = chebi_ID,
                          metlin_ID = metlin_ID, pubchem_ID = pubchem_ID, GNPS_names = gnps)
#removing "_pos" for Short_names
annotations$Short_names <- sub("_.*", "", annotations$Short_names)
#removing string starting from ";"
for (i in 2:6) {
  annotations[,i] <- gsub('([^;]*);.*', '\\1', annotations[,i])
}

#The following workflow is to utilize the web based over representation module (from enrichment analysis) to get HMDB_IDs for compounds which do not have one:
#1 Submission of compound names which do not have HMDB_ID/KEGG_ID to MetaboAnalyst web server for HMDBID conversion.
#2 Download the results and read the them into R environment.
#3 Update the HMDBID.
#4 Generate .csv file for over representation analysis
write.csv(annotations,'annotations.csv', row.names = FALSE) #Generating .csv file for compound names submission to https://www.metaboanalyst.ca/MetaboAnalyst/upload/EnrichUploadView.xhtml
#Download the result name_map.csv
MA_matched <- read.csv('path_to/name_map.csv')
match_idx <- match(annotations$Short_names, MA_matched$Query) # Find matches between Short_names and Query
annotations$HMDB_ID[!is.na(match_idx)] <- MA_matched$HMDB[match_idx[!is.na(match_idx)]] # Update HMDB_ID only where there's a match
annotations$KEGG_ID[!is.na(match_idx)] <- MA_matched$KEGG[match_idx[!is.na(match_idx)]] # Update KEGG_ID only where there's a match
annotations$pubchem_ID[!is.na(match_idx)] <- MA_matched$PubChem[match_idx[!is.na(match_idx)]] # Update Pubchem_ID only where there's a match
#Now compounds previously have no IDs are assigned one based on MetaboAnalyst database.
write.csv(annotations,'annotations_updated.csv', row.names = FALSE)

#Over representation analysis
cmpd.vec <- annotations$HMDB_ID
cmpd.vec <- unique(na.omit(cmpd.vec))
mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, cmpd.vec)
mSet<-CrossReferencing(mSet, "hmdb")
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 300, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 300, width=NA)
#For network view (i.e. .sif/.svg), refer to web version

#Quantitative Enrichment Analysis, generating the raw intensity matrix table with HMDB_ID
#The module only works for unpaired, two groups
raw_int <- read.csv('raw_intensities_matrix.csv', row.names = 1)
annotations <- read.csv('annotations_updated.csv')
matrix <- as.data.frame(metab_SE@assays@data$counts)
raw_int <- raw_int[rownames(raw_int) %in% rownames(matrix), ]
rownames(raw_int)<- metab_SE@elementMetadata$shortname
raw_int <- rbind(metab_SE@metadata$metadata$group, raw_int)
rownames(raw_int)[1] <- "Label"
replace <- c("QC","QC","QC", "Sample", "Sample", "Sample", "Sample", "Sample",
             "Sample", "Sample", "Sample", "Sample", "Sample", "Sample", "Sample",
             "Sample", "Sample", "Sample", "Sample", "Sample", "Sample")
raw_int[1,] <- replace
raw_int <- raw_int %>% 
  tibble::rownames_to_column(var = "Original_Name")
raw_int$Original_Name <- sub("_.*", "", raw_int$Original_Name)
match_idx <- match(raw_int$Original_Name, annotations$Short_names)
for (i in 1:length(raw_int$Original_Name)) {
  if (!is.na(match_idx[i])) {
    if (!is.na(annotations$HMDB_ID[match_idx[i]])) {
      raw_int$Original_Name[i] <- annotations$HMDB_ID[match_idx[i]]
    }
  }
}
duplicates <- raw_int[raw_int$Original_Name %in% raw_int$Original_Name[duplicated(raw_int$Original_Name)], ]
raw_int<- raw_int[-128, ] #Manual removal of duplicates
raw_int<- raw_int[-162, ] #Manual removal of duplicates
write.csv(raw_int, 'raw_intensity_matrix_for_qea.csv', row.names = FALSE)

mSet<-InitDataObjects("conc", "msetqea", FALSE) #For unknown reason, setting "pktable" will return 0 compound match
mSet<-Read.TextData(mSet, "raw_intensity_matrix_for_qea.csv", "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-CrossReferencing(mSet, "hmdb")
mSet<-CreateMappingResultTable(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "SumNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 300, width=NA)
mSet<-SetMetabolomeFilter(mSet, F)
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2)
mSet<-CalculateGlobalTestScore(mSet)
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "png", 300, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "png", 300, width=NA)
mSet<-SaveTransformedData(mSet)
#For network view (i.e. .sif/.svg), refer to web version