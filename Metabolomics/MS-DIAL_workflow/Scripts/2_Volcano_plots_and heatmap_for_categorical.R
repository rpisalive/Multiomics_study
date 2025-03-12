#This work flow generate volcano plots and heat maps, for categorical groups
library(RColorBrewer)

setwd('...')
metab_SE <- readRDS('.../annotated_qp.RDS') #path of the saved .rds file from script 1
source('metab_limma_categorical.R')
source('metab_limma_plot_categorical.R')
source('metab_limma_plot_heatmap.R')

# Run custom limma function on metabolites vs categorical explantory variable
metab_pos_limma_dosage <- metab_limma_categorical(metab_SE = metab_SE,
                                                      metadata_var = 'group',
                                                      metadata_condition = metab_SE@metadata$metadata$group %in% c('DIV10', 'DO50', 'DO0'), # limit data to samples from year 1
                                                      legend_metadata_string = 'Dosage_effect')

# View the volcano plot for a particular combination
metab_pos_limma_dosage$volcano_plots$`DIV10-DO0`

# Save significant values to file
#write.csv(metab_pos_limma_dosage$volcano_plots$`DIV10-DO0`,'pos_limma_dosage_DIV10DO0_significant.csv'))

# Make a list of plots for the significant breastfeeding metabolites
metab_pos_dosage_plotlist <- metab_limma_plot_all_categorical(metab_limma_cat_object = metab_pos_limma_dosage,
                                                            metab_limma_cat_comparison = 'DO0-DO50',
                                                            plot_x_lab = 'Dosage',
                                                            plot_y_lab = 'Intensity',
                                                            plot_stat_comparisons = list(c('DIV10', 'DO0'), c('DO0', 'DO50'), c('DIV10', 'DO50')),
                                                            plot_fill_name = 'Dosage',
                                                            text_size = 6.5, 
                                                            pval_text_size = 2.5)

# Arrange the top 12 metabolites in a single figure
metab_pos_bf_top12_plots <- ggarrange(plotlist = metab_pos_dosage_plotlist[1:12], nrow = 3, ncol = 4)
metab_pos_bf_top12_plots <- annotate_figure(metab_pos_bf_top12_plots, 
                                              top = text_grob('DIV10 vs. DO0 - Top 12 Differential Intensity Metabolites',
                                                              face = 'bold', size = 12))
metab_pos_bf_top12_plots

# Arrange the volcano plot and the top 12 plot next to each other
Figure_pos_bf_metab <- ggarrange(metab_pos_limma_dosage$volcano_plots$`DIV10-DO0`, metab_pos_bf_top12_plots, nrow = 1, widths = c(0.5, 0.5))
Figure_pos_bf_metab

# Define colours for extra metadata, if needed
#qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

#Original example with patient metadata
#patients <- levels(metab_stool_limma_bf_year1$input_metadata$patient)
#patient_colour_vector <- col_vector[1:length(patients)]
#names(patient_colour_vector) <- levels(metab_stool_limma_bf_year1$input_metadata$patient)
#extra_colours <- list('Patient' = patient_colour_vector)

# Run custom heatmap function
metab_pos_limma_dosage_hm <- metab_limma_plot_heatmap(metab_limma_object = metab_pos_limma_dosage, 
                                                          # metadata_to_include = c('patient'), 
                                                          # metadata_colours = extra_colours,
                                                          column_annotation_labels = 'Dosage',
                                                          categorical_colours = c('DIV10' = pal_jama(alpha = 0.6)(3)[1], # Example of setting custom colours
                                                                                  'DO50' = pal_jama(alpha = 0.6)(3)[2],
                                                                                  'DO0' = pal_jama(alpha = 0.6)(3)[3]),
                                                          heatmap_scale_name = 'Intensity',
                                                          heatmap_column_title = 'Differential Intensity Metabolites',
                                                          heatmap_row_title = 'Feature',
                                                          heatmap_show_column_names = FALSE,
                                                          heatmap_rowname_text_size = 8,
                                                          #heatmap_annotation_legend_param = list('Patient' = list(ncol = 2, by_row = TRUE)),
                                                          output_filename = ('pos_metab_dosage'))
