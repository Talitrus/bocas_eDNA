################### Differential abundances of taxa ####################

setwd("/Users/elaineshen/Documents/GitHub/bocas_eDNA/data")
library(phyloseq)
library(plyr)
library(microbiome)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)
library(data.table)

#### IMPORTING DATA ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Making RLS data into a Phyloseq file #####
  matrix = read.delim(file="RLS_fish_abundance_table.tsv",header=T, row.names = "event")
  matrix=as.matrix(matrix)
  matrix_table = otu_table(matrix,taxa_are_rows=FALSE)
  
  metadata = read.delim(file="fish_site_metadata.tsv",header=T,row.names="event")
  metadata_table = sample_data(metadata)
  
  taxa = read.csv("RLS_fish_taxonomy.csv",row.names=1)
  taxa <- taxa[,-1]
  taxmat = as.matrix(taxa)
  tax_table=tax_table(taxmat)
  
  RLS = phyloseq(matrix_table,tax_table,metadata_table)
  
  # Merge triplicate samples together, then manually change metadata
  RLS.m = merge_samples(RLS,group = "Site.Code.2") # Site.Code.2 has the habitat type followed by the site code. This helps organize the heatmap later (ex. M-ALR, Sa-ALR). 
  
  # Rename sample_variables() that got converted to #'s during the merge
  new_sample_metadata = sample_data(RLS.m)
  new_sample_metadata$Habitat <- mapvalues(new_sample_metadata$Habitat, from=c("1","2","3","4","5"), to = c("Dock","Mangrove","Reef","Sand","Seagrass"))
  new_sample_metadata$Location1 <- mapvalues(new_sample_metadata$BayRegion1, from=c("1","2"), to = c("Exposed","Sheltered"))
  sample_data(RLS.m) <- new_sample_metadata

### Importing eDNA data ###
  eDNA = readRDS(file="curated_Bocas_eDNA_phyloseq.rds")
  eDNA_metadata = read.csv(file="metadata_table_eDNA.csv",header=T,row.names = 1) # Update metadata file with one that has been manually corrected to organize heatmaps appropriately 
  sample_data(eDNA) = eDNA_metadata
  eDNA = subset_samples(eDNA, Sample=="Primary")
  eDNA.m = merge_samples(eDNA,group = "site.code4")
  # Only analyze sites with a corresponding RLS survey
  eDNA.m = prune_samples(samples=c("C-ALR",  "M-ALR",  "S-ALR",  "Sa-ALR", "D-CAR","C-CCR",  "M-CCR",  "S-CCR",  "Sa-CCR", "D-COS",  "D-FER",  "M-IPI",  "S-IPI",
                                   "M-MAR",  "S-MAR", "M-MYS",  "S-MYS",  "S-PBL",  "C-PJN",  "M-PJN",  "S-PJN",  "Sa-PJN", "C-PPR",  "M-PPR",  "S-PPR",  "Sa-PPR",
                                   "C-PST",  "D-PST",  "M-PST",  "S-PST",  "Sa-PST", "M-ROL",  "S-ROL",  "C-SCR",  "S-SCR",  "Sa-SCR","M-SCR", "S-SGL",  "D-SGN",  "M-SGN", 
                                   "S-SGN",  "M-SIS",  "S-SIS"),eDNA.m)
  # Add in Location column from RLS survey and adjust eDNA Habitat column 
  sample_data(eDNA.metazoan.m)$Location <- sample_data(RLS.m)$Location
  sample_data(eDNA.metazoan.m)$Habitat <- sample_data(RLS.m)$Habitat
  sample_data(eDNA.metazoan.m)$BayRegion1 <- sample_data(RLS.m)$Location1

##### Reef Life Survey comparisons ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### RLS - HABITAT ###

  ## RLS - Mangrove versus Seagrass ##
  RLS.m.MS = subset_samples(RLS.m, Habitat %in% c("Seagrass","Mangrove")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)
  RLS_H = phyloseq_to_deseq2(RLS.m.MS, ~Habitat) 
  RLS_H = estimateSizeFactors(RLS_H, type="poscounts")
  RLS_H = DESeq(RLS_H)

  RLS_H_vst = varianceStabilizingTransformation(RLS_H,blind=FALSE) # Variance Stabilizing Transformation uses a parametric fit for the dispersion and includes correction for size factors or normalization factors. The transformed data is on the log2 scale for large counts. This helps with visualization only and does not affect signficiance testing.

  # Determines which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05):
  RLS_H_res = results(RLS_H, contrast = c("Habitat","Mangrove","Seagrass"))
  RLS_H_res = cbind(as(RLS_H_res,"data.frame"),as(tax_table(RLS.m.MS)[rownames(RLS_H_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  RLS_H_sigtab = RLS_H_res[which(RLS_H_res$padj < alpha), ]
  RLS_H_sigtab = cbind(as(RLS_H_sigtab, "data.frame"), as(tax_table(RLS.m.MS)[rownames(RLS_H_sigtab), ], "matrix"))
  summary(RLS_H_sigtab)
  RLS_H_sigtab
  
  # Plotting the top 20 the fish taxa in order of ascending adjusted p-value (decreasing significance)
  RLS_H_select.p = order(RLS_H_res$padj)[1:20]
  # Create Habitat annotation
  Habitat = colData(RLS_H)[,c("Habitat")]
  RLS_H_df = as.data.frame(Habitat)
  rownames(RLS_H_df) <- colnames(RLS_H)
  # Manually change default RLS species code to taxa names
  RLS_taxa = as.data.frame(read.csv("RLS_fish_taxonomy.csv",row.names=1))
  RLS_taxa = RLS_taxa[order(row.names(RLS_taxa)),][,c("Species")]
  row.names(RLS_H_vst) = RLS_taxa
  
  # Set scale & aesthetics, plot
  breaksList = seq(0,20, by=1)
  
  ann_colors = list(Habitat = c(Mangrove = "magenta", Seagrass= "green"))
  pheatmap(assay(RLS_H_vst)[RLS_H_select.p,],color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList, cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=RLS_H_df, fontsize = 20, filename = "RLS_Habitat_Heatmap_top20.jpg", gaps_col=12, width=18,height=10)

  ## RLS - Sand versus Reef ## 
  RLS.m.RSa = subset_samples(RLS.m, Habitat %in% c("Sand","Reef")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)
  RLS_RSa = phyloseq_to_deseq2(RLS.m.RSa, ~Habitat) 
  RLS_RSa = estimateSizeFactors(RLS_RSa, type="poscounts")
  RLS_RSa = DESeq(RLS_RSa)
  
  RLS_RSa_vst = varianceStabilizingTransformation(RLS_RSa,blind=FALSE) # Variance Stabilizing Transformation uses a parametric fit for the dispersion and includes correction for size factors or normalizatoin factors. The transformed data is on the log2 scale for large counts.
  
  # To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  RLS_RSa_res = results(RLS_RSa, contrast = c("Habitat","Sand","Reef"))
  RLS_RSa_res = cbind(as(RLS_RSa_res,"data.frame"),as(tax_table(RLS.m.RSa)[rownames(RLS_RSa_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  RLS_RSa_sigtab = RLS_RSa_res[which(RLS_RSa_res$padj < alpha), ]
  RLS_RSa_sigtab = cbind(as(RLS_RSa_sigtab, "data.frame"), as(tax_table(RLS.m.RSa)[rownames(RLS_RSa_sigtab), ], "matrix"))
  summary(RLS_RSa_sigtab)
  RLS_RSa_sigtab
  
  # Plotting the top 20 fish taxa in order of ascending adjusted p-value (decreasing significance)
  RLS_RSa_select.p = order(RLS_RSa_res$padj)[1:20]
  
  # Create Habitat annotation
  Habitat = colData(RLS_RSa)[,c("Habitat")]
  RLS_RSa_df = as.data.frame(Habitat)
  rownames(RLS_RSa_df) <- colnames(RLS_RSa)
  # Manually change default RLS species code to taxa
  row.names(RLS_RSa_vst) = RLS_taxa
  
  breaksList = seq(0,20, by=1)
  ann_colors = list(Habitat = c(Sand = "yellow", Reef= "orange"))
  pheatmap(assay(RLS_RSa_vst)[RLS_RSa_select.p,],color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList, cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=RLS_RSa_df, fontsize = 20, filename = "RLS_Habitat_RSa_Heatmap_top20.jpg", gaps_col=6,width=12,height=10)
  
## RLS - Exposed versus Sheltered ## 
  RLS.m.L = subset_samples(RLS.m, Location1 %in% c("Exposed", "Sheltered"))
  RLS_L = phyloseq_to_deseq2(RLS.m.L, ~Location1) 
  RLS_L = estimateSizeFactors(RLS_L, type="poscounts")
  RLS_L = DESeq(RLS_L)
  
  RLS_L_vst = varianceStabilizingTransformation(RLS_L,blind=FALSE)
  
  # To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  RLS_L_res = results(RLS_L, contrast=c("Location1","Exposed", "Sheltered"))
  RLS_L_res = cbind(as(RLS_L_res,"data.frame"),as(tax_table(RLS.m)[rownames(RLS_L_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  RLS_L_sigtab = RLS_L_res[which(RLS_L_res$padj < alpha), ]
  RLS_L_sigtab = cbind(as(RLS_L_sigtab, "data.frame"), as(tax_table(RLS.m.L)[rownames(RLS_L_sigtab), ], "matrix"))
  summary(RLS_L_sigtab) 
  RLS_L_sigtab
                                                              
  # Plotting the top 20 fish taxa in order of ascending adjusted p-value (decreasing significance)
  RLS_L_select.p = order(RLS_L_res$padj)[1:20]
  
  # create Habitat annotation
  Location = colData(RLS_L)[,c("Location1")]
  RLS_L_df = as.data.frame(Location)
  rownames(RLS_L_df) <- colnames(RLS_L)
  # Manually change default RLS species code to taxa
  row.names(RLS_L_vst) = RLS_taxa
  
  # Set scale & aesthetics, plot
  breaksList = seq(0,20, by=1)
  
  pheatmap(assay(RLS_L_vst)[RLS_L_select.p,],color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList,cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_col=RLS_L_df, fontsize = 20, filename = "RLS_Location_Heatmap_top20.jpg", width=30,height=10,gaps_col=c(6,11,23,37))
  
##### eDNA comparisons ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## eDNA fish - HABITAT ## 
  
  ## eDNA fish - Mangrove versus Seagrass ##
  eDNA.fish.m = subset_taxa(eDNA.m, Class=="Actinopteri") 
  sample_data(eDNA.fish.m)$Location <- sample_data(RLS.m)$Location
  sample_data(eDNA.fish.m)$Habitat <- sample_data(RLS.m)$Habitat
  sample_data(eDNA.fish.m)$BayRegion1 <- sample_data(RLS.m)$Location1
  
  eDNA.fish.m.MG <- subset_samples(eDNA.fish.m, (Habitat %in% c("Seagrass","Mangrove"))) 
  eDNA_fish_H = phyloseq_to_deseq2(eDNA.fish.m.MG, ~Habitat) 
  eDNA_fish_H = estimateSizeFactors(eDNA_fish_H, type="poscounts")
  eDNA_fish_H = DESeq(eDNA_fish_H)
  
  eDNA_fish_H_vst = varianceStabilizingTransformation(eDNA_fish_H,blind=FALSE)
  
  ## To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_fish_H_res = results(eDNA_fish_H,contrast=c("Habitat","Mangrove","Seagrass"))
  eDNA_fish_H_res = cbind(as(eDNA_fish_H_res,"data.frame"),as(tax_table(eDNA.fish.m.MG)[rownames(eDNA_fish_H_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_fish_sigtab = eDNA_fish_H_res[which(eDNA_fish_H_res$padj < alpha), ]
  eDNA_fish_sigtab = cbind(as(eDNA_fish_sigtab, "data.frame"), as(tax_table(eDNA.fish.m.MG)[rownames(eDNA_fish_sigtab), ], "matrix"))
  summary(eDNA_fish_sigtab)
  eDNA_fish_sigtab
  
  # Plot top 20 fish taxa in order of adjusted p-value (decreasing significance)
  eDNA_fish_select.p = order(eDNA_fish_H_res$padj)[1:20]
  # create Habitat annotation
  Habitat = colData(eDNA_fish_H)[,c("Habitat")]
  eDNA_fish_df = as.data.frame(Habitat)
  rownames(eDNA_fish_df) <- colnames(eDNA_fish_H)
  # Manually change OTU ID to taxa ID (make sure both files are in alphabetical order!)
  eDNA_fish_taxa = as.data.frame(read.csv("eDNA_fish_taxonomy_table.csv",row.names=1)) # File I made with OTU ID as the rowname and lowest  taxonomic ID as the column
  ann_colors = list(Habitat = c(Mangrove = "magenta", Seagrass= "green"))
  row.names(eDNA_fish_H_vst) = eDNA_fish_taxa$ID
  breaksList = seq(0,20, by=1)
  
  pheatmap(assay(eDNA_fish_H_vst)[eDNA_fish_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList, cellwidth=25, cellheight=25, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols=FALSE,annotation_colors = ann_colors, annotation_col=eDNA_fish_df, fontsize=20, filename = "eDNA_Habitat_Heatmap_fish_20.jpg",gaps_col=12, width=18,height=10)
  
  ## eDNA fish - Sand versus Reef ##
  # I use the same variable names for metazoans for this
  eDNA.m.RSa = subset_samples(eDNA.fish.m, Habitat %in% c("Sand","Reef")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)
  eDNA_RSa = phyloseq_to_deseq2(eDNA.m.RSa, ~Habitat) 
  eDNA_RSa = estimateSizeFactors(eDNA_RSa, type="poscounts")
  eDNA_RSa = DESeq(eDNA_RSa)
  
  eDNA_RSa_vst = varianceStabilizingTransformation(eDNA_RSa,blind=FALSE) # Variance Stabilizing Transformation uses a parametric fit for the dispersion and includes correction for size factors or normalizatoin factors. The transformed data is on the log2 scale for large counts.
  
  ## To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_RSa_res = results(eDNA_RSa, contrast = c("Habitat","Sand","Reef"))
  eDNA_RSa_res = cbind(as(eDNA_RSa_res,"data.frame"),as(tax_table(eDNA.m.RSa)[rownames(eDNA_RSa_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_RSa_sigtab = eDNA_RSa_res[which(eDNA_RSa_res$padj < alpha), ]
  eDNA_RSa_sigtab = cbind(as(eDNA_RSa_sigtab, "data.frame"), as(tax_table(eDNA.m.RSa)[rownames(eDNA_RSa_sigtab), ], "matrix"))
  summary(eDNA_RSa_sigtab)
  eDNA_RSa_sigtab
  
  #Plot top 20 fish taxa in order of adjusted p-value (decreasing significance)
  eDNA_RSa_select.p = order(eDNA_RSa_res$padj)[1:20]
  # create Habitat annotation
  Habitat = colData(eDNA_RSa)[,c("Habitat")]
  eDNA_RSa_df = as.data.frame(Habitat)
  rownames(eDNA_RSa_df) <- colnames(eDNA_RSa)
  
  # Manually change OTU ID to taxa ID (make sure both files are in alphabetical order!)
  row.names(eDNA_RSa_vst) = eDNA_fish_taxa$ID
  
  breaksList = seq(0,20, by=1)
  ann_colors = list(Habitat = c(Sand = "yellow", Reef= "orange"))
  pheatmap(assay(eDNA_RSa_vst)[eDNA_RSa_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList,cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=eDNA_RSa_df, fontsize = 20, filename = "eDNA_Habitat_fish_RSa_Heatmap_20.jpg", gaps_col=6,width=12,height=10)
  
## eDNA metazoans - HABITAT ## 
  
  ## eDNA metazoans - Mangrove versus Seagrass ##
  eDNA.metazoan.m = subset_taxa(eDNA.m, Kingdom=="Metazoa") 
  eDNA.metazoan.m.MG <- subset_samples(eDNA.metazoan.m, (Habitat %in% c("Seagrass","Mangrove"))) 
  eDNA_meta_H = phyloseq_to_deseq2(eDNA.metazoan.m.MG, ~Habitat) 
  eDNA_meta_H = estimateSizeFactors(eDNA_meta_H, type="poscounts")
  eDNA_meta_H = DESeq(eDNA_meta_H)
  
  eDNA_meta_H_vst = varianceStabilizingTransformation(eDNA_meta_H,blind=FALSE)
  
  ## To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_meta_H_res = results(eDNA_meta_H,contrast=c("Habitat","Mangrove","Seagrass"))
  eDNA_meta_H_res = cbind(as(eDNA_meta_H_res,"data.frame"),as(tax_table(eDNA.metazoan.m.MG)[rownames(eDNA_meta_H_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_meta_H_sigtab = eDNA_meta_H_res[which(eDNA_meta_H_res$padj < alpha), ]
  eDNA_meta_H_sigtab = cbind(as(eDNA_meta_H_sigtab, "data.frame"), as(tax_table(eDNA.metazoan.m.MG)[rownames(eDNA_meta_H_sigtab), ], "matrix"))
  summary(eDNA_meta_H_sigtab)
  eDNA_meta_H_sigtab
  
  # Plot top 50 metazoan taxa in order of adjusted p-value (decreasing significance)
  eDNA_meta_H_select.p = order(eDNA_meta_H_res$padj)[1:50]
  # create Habitat annotation
  Habitat = colData(eDNA_meta_H)[,c("Habitat")]
  eDNA_meta_H_df = as.data.frame(Habitat)
  rownames(eDNA_meta_H_df) <- colnames(eDNA_meta_H)
  # Manually change OTU ID to taxa ID (make sure both files are in alphabetical order!)
  eDNA_meta_taxa = as.data.frame(read.csv("eDNA_taxonomy_table.csv",row.names=1)) # File I made with OTU ID as the rowname and lowest  taxonomic ID as the column
  ann_colors = list(Habitat = c(Mangrove = "magenta", Seagrass= "green"))
  row.names(eDNA_meta_H_vst) = eDNA_meta_taxa$ID
  breaksList = seq(0,20, by=1)
  
  pheatmap(assay(eDNA_meta_H_vst)[eDNA_meta_H_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList, cellwidth=25, cellheight=25, cluster_rows = FALSE, show_rownames = TRUE, cluster_cols=FALSE,annotation_colors = ann_colors, annotation_col=eDNA_meta_H_df, fontsize=20, filename = "eDNA_Habitat_Heatmap_metazoan_50.jpg", gaps_col=12,width=18,height=20)
  
  ## eDNA metazoans- Sand versus Reef ##
  eDNA.m.RSa = subset_samples(eDNA.metazoan.m, Habitat %in% c("Sand","Reef")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)
  eDNA_RSa = phyloseq_to_deseq2(eDNA.m.RSa, ~Habitat) 
  eDNA_RSa = estimateSizeFactors(eDNA_RSa, type="poscounts")
  eDNA_RSa = DESeq(eDNA_RSa)
  
  eDNA_RSa_vst = varianceStabilizingTransformation(eDNA_RSa,blind=FALSE) # Variance Stabilizing Transformation uses a parametric fit for the dispersion and includes correction for size factors or normalizatoin factors. The transformed data is on the log2 scale for large counts.
  
  ## To determine which (and the # of) metazoan taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_RSa_res = results(eDNA_RSa, contrast = c("Habitat","Sand","Reef"))
  eDNA_RSa_res = cbind(as(eDNA_RSa_res,"data.frame"),as(tax_table(eDNA.m.RSa)[rownames(eDNA_RSa_res),],"matrix"))
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_RSa_sigtab = eDNA_RSa_res[which(eDNA_RSa_res$padj < alpha), ]
  eDNA_RSa_sigtab = cbind(as(eDNA_RSa_sigtab, "data.frame"), as(tax_table(eDNA.m.RSa)[rownames(eDNA_RSa_sigtab), ], "matrix"))
  summary(eDNA_RSa_sigtab)
  eDNA_RSa_sigtab
  
  #Plot top 50 metazoan taxa in order of adjusted p-value (decreasing significance)
  eDNA_RSa_select.p = order(eDNA_RSa_res$padj)[1:50]
  # create Habitat annotation
  Habitat = colData(eDNA_RSa)[,c("Habitat")]
  eDNA_RSa_df = as.data.frame(Habitat)
  rownames(eDNA_RSa_df) <- colnames(eDNA_RSa)
  
  # Manually change OTU ID to taxa names 
  row.names(eDNA_RSa_vst) = eDNA_meta_taxa$ID

  breaksList = seq(0,20, by=1)
  ann_colors = list(Habitat = c(Sand = "yellow", Reef= "orange"))
  pheatmap(assay(eDNA_RSa_vst)[eDNA_RSa_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList,cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=eDNA_RSa_df, fontsize = 20, filename = "eDNA_Habitat_meta_RSa_Heatmap_50.jpg", gaps_col=6,width=12,height=20)

## eDNA metazoans - LOCATION ##
  
  eDNA.metazoan.m.L <- subset_samples(eDNA.metazoan.m, (BayRegion1 %in% c("Exposed","Sheltered"))) 
  eDNA_meta_L = phyloseq_to_deseq2(eDNA.metazoan.m.L, ~BayRegion1) 
  eDNA_meta_L = estimateSizeFactors(eDNA_meta_L, type="poscounts")
  eDNA_meta_L = DESeq(eDNA_meta_L)
  
  eDNA_meta_L_vst = varianceStabilizingTransformation(eDNA_meta_L,blind=FALSE)
  
  ## To determine which (and the # of) metazoan taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_meta_L_res = results(eDNA_meta_L,contrast=c("BayRegion1","Exposed","Sheltered"))
  eDNA_meta_L_res = cbind(as(eDNA_meta_L_res,"data.frame"),as(tax_table(eDNA.metazoan.m.L)[rownames(eDNA_meta_L_res),],"matrix"))
  
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_meta_L_sigtab = eDNA_meta_L_res[which(eDNA_meta_L_res$padj < alpha), ]
  eDNA_meta_L_sigtab = cbind(as(eDNA_meta_L_sigtab, "data.frame"), as(tax_table(eDNA.metazoan.m.L)[rownames(eDNA_meta_L_sigtab), ], "matrix"))
  summary(eDNA_meta_L_sigtab)
  eDNA_meta_L_sigtab

  # Plot the metazoan taxa in order of ascending p-value (decreasing significance)
  eDNA_meta_L_select.p = order(eDNA_meta_L_res$padj)[1:50]
  # create Location annotation
  Location = colData(eDNA_meta_L)[,c("BayRegion1")]
  eDNA_meta_L_df = as.data.frame(Location)
  rownames(eDNA_meta_L_df) <- colnames(eDNA_meta_L)
  # Manually change OTU ID to taxa names 
  row.names(eDNA_meta_L_vst) = eDNA_meta_taxa$ID
  
  breaksList = seq(0,20, by=1)
  pheatmap(assay(eDNA_meta_L_vst)[eDNA_meta_L_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList,cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=eDNA_meta_L_df, fontsize = 20, filename = "eDNA_Location_meta_Heatmap_50.jpg", gaps_col=c(6,11,23,37),width=25,height=20)
  
## eDNA fish - LOCATION ##
  
  eDNA.fish.m.L <- subset_samples(eDNA.fish.m, (BayRegion1 %in% c("Exposed","Sheltered"))) 
  eDNA_fish_L = phyloseq_to_deseq2(eDNA.fish.m.L, ~BayRegion1) 
  eDNA_fish_L = estimateSizeFactors(eDNA_fish_L, type="poscounts")
  eDNA_fish_L = DESeq(eDNA_fish_L)
  
  eDNA_fish_L_vst = varianceStabilizingTransformation(eDNA_fish_L,blind=FALSE)
  
  ## To determine which (and the # of) fish taxa are statsitically differnetially abundant (p<0.05), use:
  eDNA_fish_L_res = results(eDNA_fish_L,contrast=c("BayRegion1","Exposed","Sheltered"))
  eDNA_fish_L_res = cbind(as(eDNA_fish_L_res,"data.frame"),as(tax_table(eDNA.fish.m.L)[rownames(eDNA_fish_L_res),],"matrix"))
  
  # (Optional) A list of the statistically significant species (defined by your alpha value)
  alpha = 0.05
  eDNA_fish_L_sigtab = eDNA_fish_L_res[which(eDNA_fish_L_res$padj < alpha), ]
  eDNA_fish_L_sigtab = cbind(as(eDNA_fish_L_sigtab, "data.frame"), as(tax_table(eDNA.fish.m.L)[rownames(eDNA_fish_L_sigtab), ], "matrix"))
  summary(eDNA_fish_L_sigtab)
  eDNA_fish_L_sigtab
  
  # Plot the top 20 fish taxa in order of ascending p-value (decreasing significance)
  eDNA_fish_L_select.p = order(eDNA_fish_L_res$padj)[1:20]
  # create Location annotation
  Location = colData(eDNA_fish_L)[,c("BayRegion1")]
  eDNA_fish_L_df = as.data.frame(Location)
  rownames(eDNA_fish_L_df) <- colnames(eDNA_fish_L)
  # Manually change OTU ID to taxa names 
  row.names(eDNA_fish_L_vst) = eDNA_fish_taxa$ID
  
  breaksList = seq(0,20, by=1)
  pheatmap(assay(eDNA_fish_L_vst)[eDNA_fish_L_select.p,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),breaks=breaksList,cellwidth = 25, cellheight=25, cluster_rows = FALSE,show_rownames = TRUE, cluster_cols=FALSE, annotation_colors = ann_colors, annotation_col=eDNA_fish_L_df, fontsize = 20, filename = "eDNA_Location_fish_Heatmap_50.jpg", gaps_col=c(6,11,23,37),width=24,height=15)
  
