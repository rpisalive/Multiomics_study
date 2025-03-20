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
#4 Generate .csv file for over representation analysis. When submitting online, make sure to use both HMDBID and KEGGID for maximum coverage!
#Pure R version is NOT recommeneded due to the usage of older vertison of HMDB.
write.csv(annotations,'annotations.csv', row.names = FALSE) #Generating .csv file for compound names submission.
converted <- read.csv('path_to/name_map.csv') #Read the results into the environment.
match_idx <- match(annotations$Short_names, converted$Query) # Find matches between Short_names and Query
annotations$HMDB_ID[!is.na(match_idx)] <- converted$HMDB[match_idx[!is.na(match_idx)]] # Update HMDB_ID only where there's a match
annotations$KEGG_ID[!is.na(match_idx)] <- converted$KEGG[match_idx[!is.na(match_idx)]] # Update KEGG_ID only where there's a match
annotations$pubchem_ID[!is.na(match_idx)] <- converted$PubChem[match_idx[!is.na(match_idx)]] # Update Pubchem_ID only where there's a match
#Now compounds previously have no IDs are assigned one based on MetaboAnalyst database.
write.csv(annotations,'annotations_updated.csv', row.names = FALSE) #Generating .csv file for over representation analysis.

#For HMDB_ID not identified by web version MetaboAnalyst
search_result <- read.csv('path_to/name_map.csv')
not_reg <- search_result[is.na(search_result$Match), ]
match_idx <- match(not_reg$Query, annotations$HMDB_ID)
not_reg$KEGG <- annotations$KEGG_ID[match_idx]
not_reg$PubChem <- annotations$pubchem_ID[match_idx]
write.csv(not_reg, 'query_for_unrecognized.csv')
#It seems like if HMDB_ID has not matches, neither would KEGG_ID or PubChem_ID.

#MetaboAnalyst database matching, based on names, HMDB ID, KEGG ID and GNPS names
#R code version of the previous section.
#However, web version is recommended because the R code is applies an older version of HMDB.
mSet_list <- list()
for (i in 1: ncol(annotations)) {
  mSet<-InitDataObjects("pktable", "msetora", TRUE)
  mSet<-Setup.MapData(mSet, annotations[,i])
  if (colnames(annotations)[i] == "Short_names" | colnames(annotations)[i] == "GNPS_names") {
    mSet<-CrossReferencing(mSet, "name", chebi = T, metlin = T, lipid = T)
  } else if (colnames(annotations)[i] == "HMDB_ID") {
    mSet<-CrossReferencing(mSet, "hmdb", chebi = T, metlin = T, lipid = T)
  } else if (colnames(annotations)[i] == "KEGG_ID") {
    mSet<-CrossReferencing(mSet, "kegg", chebi = T, metlin = T, lipid = T)
  } else if (colnames(annotations)[i] == "chebi_ID") {
    mSet<-CrossReferencing(mSet, "chebi", chebi = T, metlin = T, lipid = T)
  } else if (colnames(annotations)[i] == "metlin_ID") {
    mSet<-CrossReferencing(mSet, "metlin", chebi = T, metlin = T, lipid = T)
  }else if (colnames(annotations)[i] == "pubchem_ID") {
    mSet<-CrossReferencing(mSet, "pubchem", chebi = T, metlin = T, lipid = T)
  } else {
    stop("Unknown annotation source! Please check again.")
  }
  mSet_list[[i]] <- mSet
}
#Merge all cross referencing results into the same mSet object
for (i in 1:length(mSet_list)) {
  if (class(mSet_list[[i]]) != "list") {
    mSet_list[[i]] <- NULL
  }
}
mSet <- mSet_list[[1]]  # Copy structure from first mSet object
vector_length <- length(mSet_list[[1]]$name.map$query.vec)
mSet$name.map$match.state <- rep(0, vector_length)
mSet$name.map$hit.inx <- rep(NA, vector_length)
mSet$name.map$hit.values <- rep(NA, vector_length)
# Iterate over each index position in the vectors
for (j in 1:vector_length) {
  match_found <- FALSE  # Flag to check if we find a match
  for (i in 1:length(mSet_list)) {
    if (mSet_list[[i]]$name.map$match.state[j] == 1) {
      # Update new_mSet with the values from the first found match
      mSet$name.map$match.state[j] <- 1
      mSet$name.map$hit.inx[j] <- mSet_list[[i]]$name.map$hit.inx[j]
      mSet$name.map$hit.values[j] <- mSet_list[[i]]$name.map$hit.values[j]
      match_found <- TRUE
      break  # Stop checking once a match is found
    }
  }
  # If no match was found in any object, keep default NA values
  if (!match_found) {
  mSet$name.map$hit.inx[j] <- NA
    mSet$name.map$hit.values[j] <- NA
    mSet$name.map$match.state[j] <- 0
  }
}

#Quantitative Enrichment Analysis, generating the .csv file for web submission, only works for unpaired, two groups
raw_int <- read.csv('raw_intensities_matrix.csv', row.names = 1)
metab_SE <- readRDS('annotated_qp.RDS')
annotations <- read.csv('annotations_updated.csv')
matrix <- as.data.frame(metab_SE@assays@data$counts)
raw_int <- raw_int[rownames(raw_int) %in% rownames(matrix), ]
rownames(raw_int)<- metab_SE@elementMetadata$shortname
raw_int <- rbind(metab_SE@metadata$metadata$group, raw_int)
rownames(raw_int)[1] <- "Label"
#Please modify the code to define labels according to your experiment, refer to https://dev.metaboanalyst.ca/docs/Format.xhtml
replace <- c("Pre","Pre","Pre", "Dosed","Dosed","Dosed","Pre","Pre","Pre",
             "Pre","Pre","Pre","Dosed","Dosed","Dosed",
             "Pre","Pre","Pre","Pre","Pre","Pre") #Edut group after submission to web!
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