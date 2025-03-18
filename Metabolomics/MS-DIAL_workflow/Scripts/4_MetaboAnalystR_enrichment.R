library(MetaboAnalystR)
#library(httr)

setwd('C:/Users/49152/Downloads/Brain_trauma')
metab_SE <- readRDS('annotated_qp.RDS')
# Get HMDB and KEGG annotations
shortnames <- rowData(metab_SE)$shortname
hmdb_accession <- rowData(metab_SE)$HMDB_accession
kegg_annotations <- rowData(metab_SE)$KEGG
gnps <- rowData(metab_SE)$compound_name_gnps
annotations <- data.frame(Short_names = shortnames, HMDB_ID = hmdb_accession, KEGG_ID = kegg_annotations, GNPS_names = gnps)
#removing "_pos" for Short_names
annotations$Short_names <- sub("_.*", "", annotations$Short_names)
#removing string starting from ";"
annotations$HMDB_ID <- gsub('([^;]*);.*', '\\1', annotations$HMDB_ID)
annotations$KEGG_ID <- gsub('([^;]*);.*', '\\1', annotations$KEGG_ID)
annotations$GNPS_names <- gsub('([^;]*);.*', '\\1', annotations$GNPS_names)

#The following workflow is to utilize the web based over representation module (from enrichment analysis) to get HMDB_IDs for compounds which do not have one:
#1 Submission of compound names which do not have HMDB_ID/KEGG_ID to MetaboAnalyst web server for HMDBID conversion.
#2 Download the results and read the them into R environment.
#3 Update the HMDBID.
#4 Generate .csv file for over representation analysis. When submitting online, make sure to use both HMDBID and KEGGID for maximum coverage!
#Pure R version is NOT recommeneded due to the usage of older vertison of HMDB.
write.csv(annotations,'annotations.csv', row.names = FALSE) #Generating .csv file for compound names submission.
converted <- read.csv('name_map.csv') #Read the results into the environment.
match_idx <- match(annotations$Short_names, converted$Query) # Find matches between Short_names and Query
annotations$HMDB_ID[!is.na(match_idx)] <- converted$HMDB[match_idx[!is.na(match_idx)]] # Update HMDB_ID only where there's a match
annotations$KEGG_ID[!is.na(match_idx)] <- converted$KEGG[match_idx[!is.na(match_idx)]] # Update HMDB_ID only where there's a match
#Now compound names have IDs assigned based on MetaboAnalyst data base.
#Next, 
write.csv(annotations,'annotations_updated.csv', row.names = FALSE) #Generating .csv file for over representation analysis.

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
  } else {
    stope("Unknown annotation source! Please check again.")
  }
  mSet_list[[i]] <- mSet
}
#Merge all cross referencing results into the same mSet object
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

#The following workflow is for Compound name mapping function, currently not working due to server internal error
#name.vec <- annotations$HMDB_ID
#toSend = list(queryList = name.vec, inputType = "hmdb")
#call <- "https://rest.xialab.ca/api/mapcompounds"
#query_results <- httr::POST(call, body = toSend, encode = "json")
#query_results$status_code==500
#query_results_text <- content(query_results)
#query_results_json <- rjson::fromJSON(query_results_text, flatten = TRUE)
#query_results_table <- t(rbind.data.frame(query_results_json))
#rownames(query_results_table) <- query_results_table[,1]

#Quantitative Enrichment Analysis, generating the .csv file for web submission, only works for unpaired, two groups
#It seems that compounds which have KEGG_ID can be addressed by HMDB_ID already, so a search of HMDB_ID on web server should be sufficient, please always check though.
annotations <- read.csv('annotations_updated.csv')
raw_int_matrix <- read.csv('C:/Users/49152/Downloads/Brain_trauma/MetaboAnalyst_stat_output/meta_stat.csv') #meta_stat.csv generated for statistical analysis workflow.
raw_int_matrix$X <- sub("_.*", "", raw_int_matrix$X)
match_idx <- match(raw_int_matrix$X, annotations$Short_names)
for (i in 1:length(raw_int_matrix$X)) {
  if (!is.na(match_idx[i])) {
    if (!is.na(annotations$HMDB_ID[match_idx[i]])) {
      raw_int_matrix$X[i] <- annotations$HMDB_ID[match_idx[i]]
    } else if (!is.na(annotations$KEGG_ID[match_idx[i]])) {
      raw_int_matrix$X[i] <- annotations$KEGG_ID[match_idx[i]]
    }
  }
}
duplicates <- raw_int_matrix[raw_int_matrix$X %in% raw_int_matrix$X[duplicated(raw_int_matrix$X)], ]
raw_int_matrix<- raw_int_matrix[-128, ] #Manual removal of duplicates
raw_int_matrix<- raw_int_matrix[-162, ] #Manual removal of duplicates
replace <- c("Label","Dosed","Dosed","Dosed","Pre","Pre","Pre","Pre","Pre","Pre","Dosed","Dosed","Dosed","Pre","Pre","Pre","Pre","Pre","Pre")
raw_int_matrix[1,] <- replace
write.csv(raw_int_matrix, 'raw_intensity_matrix_for_qea.csv', row.names = FALSE)
