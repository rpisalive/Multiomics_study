gnps_format_names <- function(gnps_pos_df, gnps_neg_df) {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse', 'SummarizedExperiment', 'S4Vectors', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  gnps_remove_tags <- function(vector) {
    vector <- gsub('Spectral Match to (.*) from NIST14', '\\1', vector) # fix "spectal match" tags
    vector <- gsub('ReSpect:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "ReSpect" tags
    vector <- gsub('Massbank:[a-zA-Z0-9]* ([^|]*).*', '\\1', vector) # fix "Massbank" tags
    vector <- gsub('MoNA:\\d* (.*)', '\\1', vector) # fix "MoNA" tags
    vector <- gsub('(.*) - \\d{1,3}.\\d{1,3} eV', '\\1', vector) # fix "eV" tags
  }
  
  if (is.null(gnps_pos_df) & is.null(gnps_neg_df) ) {
    stop('Both gnps_pos_df & gnps_neg_df are NULL!')
  } else if (is.null(gnps_pos_df) | is.null(gnps_neg_df)) {
    if (is.null(gnps_neg_df)) {
      gnps_neg_df <- gnps_pos_df
    } else {
      gnps_pos_df <- gnps_neg_df
    }
  }
  
  # Reorder the columns to put #SCAN# first
  scan_colname_pos <- colnames(gnps_pos_df)[str_detect(colnames(gnps_pos_df), '#S')]
  scan_colname_neg <- colnames(gnps_neg_df)[str_detect(colnames(gnps_neg_df), '#S')]
  
  gnps_pos_df <- gnps_pos_df %>% dplyr::select(all_of(scan_colname_pos), everything())
  gnps_neg_df <- gnps_neg_df %>% dplyr::select(all_of(scan_colname_neg), everything())
  
  # Change first column name from '#SCAN#'
  colnames(gnps_pos_df)[1] <- 'Alignment.ID'
  colnames(gnps_neg_df)[1] <- 'Alignment.ID'
  
  # Copy Compound Name column
  gnps_pos_df$compound_name_gnps <- gnps_pos_df$Compound_Name
  gnps_neg_df$compound_name_gnps <- gnps_neg_df$Compound_Name
  
  # Remove the extraneous text around metabolite feature names
  gnps_pos_df$compound_name_gnps <- gnps_remove_tags(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- gnps_remove_tags(gnps_neg_df$compound_name_gnps)
  
  # Switch to title case
  gnps_pos_df$compound_name_gnps <- stringr::str_to_title(gnps_pos_df$compound_name_gnps)
  gnps_neg_df$compound_name_gnps <- stringr::str_to_title(gnps_neg_df$compound_name_gnps)
  
  # Create vector for new rownames to match pmp rownames
  new_rownames <- c(paste0(gnps_pos_df$Alignment.ID, '_pos'), paste0(gnps_neg_df$Alignment.ID, '_neg'))
  
  if (!is.null(gnps_pos_df) & !is.null(gnps_neg_df)) {
    # rbind the data.frames
    gnps_df <- rbind(gnps_pos_df, gnps_neg_df) %>%
      dplyr::select(Alignment.ID, compound_name_gnps, everything()) %>%
      mutate(alignment_ionisation = new_rownames)
    rownames(gnps_df) <- new_rownames
  } else {
      gnps_df <- gnps_pos_df
    }
  gnps_df
}