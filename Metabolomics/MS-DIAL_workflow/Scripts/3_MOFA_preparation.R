setwd('...')
metab_SE <- readRDS('annotated_qp.RDS')

#Genration of intensity matrix (normalized) for MOFA
matrix <- metab_SE@assays@data$counts
rownames(matrix) <- metab_SE@elementMetadata$shortname
#write.csv(matrix, "test.csv", row.names = TRUE)
