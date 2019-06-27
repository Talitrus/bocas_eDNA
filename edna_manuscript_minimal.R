#!/usr/bin/env Rscript

# Flags ----------------------

use_plotly_cloud <- FALSE # set this to TRUE to enable the use of `api_create` to make/update Plot.ly cloud figures.
plotly_sharing_setting <- "secret" # Plot.ly cloud privacy setting for figures generated, only applicable if `use_plotly_cloud` == TRUE.
use_h2o <- TRUE # Use `h2o` instead of the R `randomForest` package?
use_h2o_grid_search <- FALSE # Use a full grid search instead of a single model for h2o, only applicable if `use_h2o` == TRUE.
set.seed(7270) # Setting a seed for consistency.
h2o_mem_allotted <- "2g" # Set the maximum memory to be used by h2o. Better larger than smaller, if possible.

# Load libraries ------------------------

library(phyloseq)
library(vegan)
library(tidyverse)
library(plotly)
library(lubridate)
library(genefilter) # for the kOverA function
library(taxize)
if (use_h2o == TRUE) {
  library(h2o)
} else {
  library(randomForest)
}
library(DESeq2)
library(extrafont) # I make use of the free font Open Sans, so you'll need that to have the code run as-is. Feel free to edit out the custom font for your own purposes.

# Functions --------------
deseq_wrapper <- function(ps, formula = ~1, sfType = "poscounts", fitType = "parametric", test = "Wald") { # wrapper, phyloseq object as first arg
  require(phyloseq)
  pruned_ps <- prune_samples(setdiff(sample_names(ps), names(which(rowSums(otu_table(ps)) == 0))), ps) # Remove samples with zero reads
  DS_obj <- phyloseq_to_deseq2(pruned_ps, formula)
  DS_obj <- DESeq(DS_obj, test = test, fitType = fitType, sfType = sfType)
  return(DS_obj)
}

res_to_tsv <- function(x, filename, alpha_threshold = 1, ...) { #Convert DESeq results to TSV
  require(tibble)
  require(dplyr)
  require(readr)
  export_tb <- tibble::as_tibble(x@listData)
  export_tb$gene <- x@rownames
  export_tb <- export_tb %>%
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    dplyr::arrange(padj)
  if (alpha_threshold > 0) {
    export_tb <- export_tb %>%
      dplyr::filter(padj <= alpha_threshold)
  }
  readr::write_tsv(export_tb, path = filename, ...)
}

gm_mean = function(x, na.rm=TRUE){ #geometric mean function that can deal with zeroes.
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


vegan_otu <- function(physeq) { #convert phyloseq OTU table into vegan OTU matrix
  require(phyloseq)
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

filter_taxa_to_other <- function(physeq, filterFunc, merged_taxon_name = "Other taxa") { #returns phyloseq object
  require(phyloseq)
  taxa_to_merge <- (!filter_taxa(physeq, filterFunc, FALSE)) # named logical vector where samples NOT meeting the filter threshold are TRUE
  taxa_to_merge_names <- which(taxa_to_merge == TRUE) %>% names() #get the names
  new_physeq <- merge_taxa(physeq, taxa_to_merge_names, archetype = 1)
  tax_table(new_physeq)[taxa_to_merge_names[1],] <- replicate(ncol(tax_table(new_physeq)), merged_taxon_name) # merges with first taxon, so we can replace the lineage information for the first taxon to whatever we want displayed
  return(new_physeq)
}

plotly_permanova <- function(x, ...) { # Make PERMANOVA output into a plot.ly bar chart.
  require(plotly)
  require(tibble)
  if(is(x,"anova")) {
    x.tb <- x %>%
      rownames_to_column(var = "term") %>%
      as_tibble() %>%
      mutate(text = paste0("Df = ", Df, "\n F = ", round(F, digits = 2), "\n P = ", `Pr(>F)`)) %>%
      filter(term != "Total")
    
    perm_plot <- plot_ly(x.tb, y = ~term, x =~ R2, text = ~ text, textposition = "auto", type = "bar", ...) %>%
      layout(yaxis = list(title = "Term"), xaxis = list(title = "R<sup>2</sup>"))
      return(perm_plot)
  } else {
    stop("Function plotly_permanova() requires an 'anova' class object as input.")
  }
  
}

lowest_taxonomy_from_table <- function(x, tax_table) { # x should be a sha1 name or taxon name, tax_table is the phyloseq taxonomy table
  return(as.character(tax[x,ncol(tax)-sum(is.na(tax[x,]))]))
}

plotUpDownSigTaxa_eco <- function(results, rld, metadata_cols = "Ecosystem", ntaxa = 20, title, type = "sig", alpha = 1, ...) { # derived from https://www.biostars.org/p/178748/#178873, colNums normally is used to specify which samples to include, rld is the (transformed) data, results is a DESeq2 result object, ngenes is how many genes to include, metadata cols is a vector of strings matching the names of the columns to display from the metadata.

  require(SummarizedExperiment)
  require(tidyr)
  require(tibble)
  
  # make the lists
  if (alpha < 1) {
    print("Running alpha threshold.")
    results <- results[!is.na(results$padj),]
    results <- results[results$padj < alpha,]
  }
  if (type == "sig") {
    select_taxa <- rownames(head(results[ order( results$padj ), ], n=ntaxa))
    results_subset <- results[select_taxa,]
    #select_taxa <- rownames(results[order(results$log2FoldChange),])
  } else if (type == "down") {
    select_taxa <- rownames(head(results[ order( results$log2FoldChange ), ], n=ntaxa))
    results_subset <- results[select_taxa,]
    select_taxa <- rownames(results[order(results$log2FoldChange),])
  } else if (type == "up") {
    select_taxa <- rownames(head(results[ order( -results$log2FoldChange ), ], n=ntaxa))
    results_subset <- results[select_taxa,]
    select_taxa <- rownames(results[order(results$log2FoldChange),])
  } else if (type == "magnitude") {
    select_taxa <- rownames(head(results[ order( -abs(results$log2FoldChange )), ], n=ntaxa))
    results_subset <- results[select_taxa,]
    select_taxa <- rownames(results_subset[order(results_subset$log2FoldChange),])
  } else if (type == "none") {
    select_taxa <- rownames(results)
    results_subset <- results
  } else {
    stop('Type argument not recognized. Type must be one of: "sig", "down", "up", or "magnitude".')
  }
  
  # this gives us the rows we want
  rows <- match(select_taxa, row.names(rld))
  mat <- cbind(df,t(assay(rld)[rows,]))
  #mat <- mat - rowMeans(mat) # BN: Do we need this line? Try with and without?
  mat_tb <- gather(as_tibble(mat, rownames = "Site", .name_repair = "unique"), key = Taxon, value = Abundance, -c(1:2)) %>% arrange(Ecosystem)
  mat_tb$Abundance[mat_tb$Abundance < 0] <-  0 # Values < 0 set to zero.
  print("Matrix 1 generated.")
  mat_tb2 <- mat_tb
  mat_tb2$Site <- factor(mat_tb$Site, levels = unique(mat_tb$Site)) #
  mat_tb2$Taxon <- factor(mat_tb2$Taxon, levels = select_taxa)
  #mat_tb2$lowest_taxon <- sapply(as.character(mat_tb2$Taxon), lowest_taxonomy_from_table, tax_table = tax)
  # the labels are hard coded at the moment :(
  #df <- as.data.frame(colData(rld)[metadata_cols])
  heatmap1 <- ggplot(data = mat_tb2, mapping = aes(x = Site, y = Taxon, fill = Abundance)) + 
    geom_tile() +
    viridis::scale_fill_viridis() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(heatmap1)
}


# Set file locations -----------------------

# Locations to read from:
RLS_survey_location <- "data/180116_Bocas_RLS_fish-biomass_full_M1_SizeCorrected.csv" # RLS fish diver survey file location
eDNA_phyloseq_location <- "data/curated_Bocas_eDNA_phyloseq.rds" #DADA2-processed, LULU-curated OTUs in phyloseq format
fish_site_metadata_location <- "data/fish_site_metadata.tsv" #Metadata for RLS survey sites
fish_taxonomy_location <- "data/RLS_fish_taxonomy.csv" #Taxonomic lineage info file for RLS fish taxa
unprocessed_bocasDB_file <- "data/bocas_db_metazoa_20190307.txt" # Bocas species DB (metazoans only) as manually scraped from the STRI website.
GISD_data_location <- "data/GISD_04112019.csv" # Global invasive species database search for Caribbean invasives.

# Locations to write to:
metazoan_taxa_importance_file_location <- "output/randomforest_predictor_importance.txt" # Random forest feature importance file location

processed_bocasDB_file <- "output/bocasDB_metazoan_processed_20190307.tsv" # Reformatted Bocas species database, with additional information added on by this script.
#if processed_bocas_DB_files already exists, it will not be overwritten. Instead, it will be loaded and used to save time.

eDNA_spp_file <- "output/eDNA_metazoan_species.tsv" # species-level identifications from eDNA. Loaded and not overwritten if already exists.


# Load phyloseq --------------------
curated_bocas_edna.ps <- readRDS(eDNA_phyloseq_location) # DADA2, LULU-processed phyloseq object
tax <- tax_table(curated_bocas_edna.ps)
sample_data_renamed <- sample_data(curated_bocas_edna.ps) %>% # we create this file to manually preserve some of the metadata to use after split samples are merged in the next chunk.
  as.data.frame() %>%
  remove_rownames() %>%
  filter(!duplicated(site_code3)) %>%
  column_to_rownames(var = "site_code3") # Rownames are set to the "site_code3", e.g. CCR-SA, PST-RB since this is what the new sample names become after merging


# Transform data -----------------------------------
meta_curated_bocas_edna.ps.tr <- curated_bocas_edna.ps %>%
  subset_samples(Fraction == "eDNA" & Sample == "Primary") %>% #filter out control, negative and mock samples # if desired, we can also use  `& Ecosystem != "Dock"` to remove docks
  merge_samples("site_code3", fun = mode) %>% # Doesn't matter what the `fun` argument is set to, I think, we overwrite it in the next step
  subset_taxa(Kingdom == "Metazoa") %>% # Metazoan only
  filter_taxa(function (x) sum(x) > 0, TRUE) %>% # Remove taxa that are no longer represented because of subsetting
  transform_sample_counts(function(x) x / sum(x)) # Transform to relative abundances

sample_data(meta_curated_bocas_edna.ps.tr) <- sample_data(sample_data_renamed)

eDNA_meta.ps <- curated_bocas_edna.ps %>% # same as above but not transformed to relative abundances
  subset_samples(Fraction == "eDNA" & Sample == "Primary") %>% #filter out control, negative and mock samples # if desired, we can also use  `& Ecosystem != "Dock"` to remove docks
  merge_samples("site_code3", fun = mode) %>% # Doesn't matter what the `fun` argument is set to, I think, we overwrite it in the next step
  subset_taxa(Kingdom == "Metazoa") %>% # Metazoan only
  filter_taxa(function (x) sum(x) > 0, TRUE)# Remove taxa that are no longer represented because of subsetting

sample_data(eDNA_meta.ps) <- sample_data(sample_data_renamed)

# Load fish data -------------------------------------

fish_data <- read_csv(RLS_survey_location) %>%
  filter(Diver != 'Daniella Heflin') %>% # We had observations from a third diver for some sites—removed for consistency
  select(Line_ID, Diver, Buddy, Site.No., event, Site.Name, Ecosystem, Latitude, Longitude, Date, vis, Direction, Time, Depth, Method, Block, Code, Species, Common.name, `Total Abundance`, `Total Biomass`) %>%
  mutate(Code = str_to_lower(Code)) %>% # Code is an abbrev. for species, lowercased for consistency
  unique() # One row is an exact duplicate, likely an error. Removing duplicates
fish_data$Species[fish_data$Species == "Clupeid spp."] <- "Clupeidae" # manual correction

if(file.exists(fish_site_metadata_location)) { #If a metadata file exists, read it in
  # Currently, this is used because I manually added a column called BayRegion to the metadata file, which corresponds to "BayRegion" in the eDNA
  # At some point, Location should be renamed to BayRegion for consistency
  print("Fish site metadata exists. Loading in fish site metadata.")
  fish_site_metadata_sub <- read.delim(fish_site_metadata_location) %>% # Load in a new metadata file
    mutate(Site = toupper(str_remove(Site.Name, " (Sand|Reef|Mangroves?|Dock2?|Seagrass)")) %>%
             str_replace("STRI POINT", "STRI") %>%
             str_replace("SAIGON BAY", "SAIGON")) %>%
    column_to_rownames(var = "event")
} else { # Otherwise, make a metadata file from the RLS survey data file.
  fish_site_metadata <- fish_data %>%
    unite(Date, Time, col = "datetime", sep = " ") %>%
    mutate(
      datetime = mdy_hms(datetime, tz = "America/Panama"), 
      Site = toupper(str_remove(Site.Name, " (Sand|Reef|Mangroves?|Dock2?|Seagrass)")) %>% # Cleaning to extract site names
        str_replace("STRI POINT", "STRI") %>% # for consistency
        str_replace("SAIGON BAY", "SAIGON")
    ) %>%
    select(event, Site.Name,Ecosystem, Latitude, Longitude, Depth, Method, Block, datetime, Diver) %>%
    unique()
  fish_site_metadata_sub <- column_to_rownames(as.data.frame(fish_site_metadata), var = "event")
}

fish_code_dictionary <- fish_data %>% #Used for looking up which RLS code corresponds to which scientific name, for manual checking only—not used further in this code
  select(Code, Species) %>%
  unique() %>%
  filter((Code != 'hst') | (Species != "Haemulon parra"))

# Converting RLS divery survey data entry format to a species abundance table
fish_spread <- spread(fish_data %>% select(event, Code, `Total Abundance`), key = 'Code', value = 'Total Abundance', fill = 0)
fish_table <- as.data.frame(fish_spread[,-1])
rownames(fish_table) <- fish_spread$event


# Load fish taxonomy table

fish_tax <- read_csv(file = fish_taxonomy_location) %>%
  arrange(Code) %>%
  select(-Species)
fish_tax <- as.matrix(column_to_rownames(fish_tax, var = "Code"))


# how many fish IDs to spp level out of total observations?

sum(!is.na(fish_tax[match(fish_data$Code, rownames(fish_tax)),"species"])) # number of species-level ID events
nrow(fish_data) # total number of observation events.

fish_ps <- phyloseq(otu_table(as.matrix(fish_table), taxa_are_rows = FALSE), tax_table(fish_tax), sample_data(fish_site_metadata_sub))

# Taxonomic composition -------------------------------------------

# Superkingdom-level taxonomic statistics
skglom_curated_bocas_edna.ps.tr <- curated_bocas_edna.ps %>% # Aggregate data into phylum-level.
  subset_samples(Fraction == "eDNA") %>% #filter out control, negative and mock samples, NOTE: currently includes non-primary samples
  tax_glom(taxrank = "Superkingdom")

# How many reads from each superkingdom?
(taxa_sums(skglom_curated_bocas_edna.ps.tr))/sum(sample_sums(curated_bocas_edna.ps))
tax[names(taxa_sums(skglom_curated_bocas_edna.ps.tr)),]

# Kingdom-level taxonomic statistics
kingglom_curated_bocas_edna.ps.tr <- curated_bocas_edna.ps %>% # Aggregate data into phylum-level.
  subset_samples(Fraction == "eDNA") %>% #filter out control, negative and mock samples, NOTE: currently includes non-primary samples
  tax_glom(taxrank = "Kingdom")

# How many reads from each superkingdom?
(taxa_sums(kingglom_curated_bocas_edna.ps.tr))/sum(sample_sums(curated_bocas_edna.ps))
tax[names(taxa_sums(kingglom_curated_bocas_edna.ps.tr)),]

#How many metazoan reads and OTUs?
# How many Metazoan OTUs?
meta_counts.ps <- curated_bocas_edna.ps %>%
  subset_samples(Fraction == "eDNA" & Sample == "Primary") %>% # Doesn't matter what the `fun` argument is set to, I think, we overwrite it in the next step
  subset_taxa(Kingdom == "Metazoa")

sum(sample_sums(meta_counts.ps))
ntaxa(meta_counts.ps) 

# How many reads identified down to species?
# See end of script.



# Phylum-level taxonomic composition plot
phyglom_curated_bocas_edna.ps.tr <- curated_bocas_edna.ps %>% # Aggregate data into phylum-level.
  subset_samples(Fraction == "eDNA") %>% #filter out control, negative and mock samples, NOTE: currently includes non-primary samples
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  filter_taxa_to_other(function(x) max(x) > 0.03) # All taxa that are never > 0.03 rel. abun. merged into "Other taxa"

edna_phylum <- plot_ly( data = psmelt(phyglom_curated_bocas_edna.ps.tr), x = ~Abundance, y = ~Sample, color =~ Phylum, type = 'bar', colors = "Paired") %>%
  layout(barmode = "stack",
         bargap = 0,
         legend = list( traceorder = "normal"),
         xaxis = list(title = "Relative Abundance",
                      range = c(0,1),
                      showgrid = FALSE,
                      ticks = "outside"
                      ), 
         yaxis = list(title = "Samples",
                      categoryorder = "category ascending",
                      showticklabels = FALSE,
                      showgrid = FALSE
                      )
         )
edna_phylum

if (use_plotly_cloud == TRUE) {
  api_create(edna_phylum, filename = "bocas/edna/tax_comp/edna_phylum", sharing = plotly_sharing_setting)
}


# Random Forest habitat classification --------------------

# Pre-processing the input data
curated_edna_meta.tr.f <- meta_curated_bocas_edna.ps.tr %>%
  filter_taxa(filterfun(kOverA(5,5e-5)), TRUE) #filter for taxa with approximately at least 5 reads per sample in at least 5 samples)
selected_env <- sample_data(curated_edna_meta.tr.f) %>% # column containing the class information (Ecosystem, in this case)
  select(Ecosystem)
meta_comb <- cbind(otu_table(curated_edna_meta.tr.f),selected_env) #combined dataframe with OTUs plus the Ecosystem column.

if (use_h2o == FALSE) {
  meta_rf <- randomForest(meta_comb[,1:(ncol(meta_comb)-1)], meta_comb[,ncol(meta_comb)], ntree = 2000, importance = TRUE, proximity = TRUE, mtry = 200) # Random forest algorithm. Tweak settings here if necessary—if not using h2o.
  raw_importance <- importance(meta_rf) %>% # Feature importance
    as.data.frame() %>%
    rownames_to_column(var = "SHA1") %>% 
    as_tibble() 
  
  meta_imp_table <- raw_importance %>% # Appending taxonomy information
    mutate(
      Superkingdom = as.character(tax[raw_importance$SHA1,"Superkingdom"]),
      Kingdom = as.character(tax[raw_importance$SHA1,"Kingdom"]), 
      Phylum = as.character(tax[raw_importance$SHA1,"Phylum"]),
      Class = as.character(tax[raw_importance$SHA1,"Class"]),
      Order = as.character(tax[raw_importance$SHA1,"Order"]),
      Family = as.character(tax[raw_importance$SHA1,"Family"]),
      Genus = as.character(tax[raw_importance$SHA1,"Genus"]),
      Species = as.character(tax[raw_importance$SHA1,"species"])
    ) %>%
    arrange(desc(MeanDecreaseAccuracy)) # Sort features by decreasing Mean Decrease in Accuracy
  
  write_tsv(meta_imp_table, path = metazoan_taxa_importance_file_location) # Write variable importances to a file.
}




# H2o distributed random forest ----------------

if (use_h2o == TRUE) { # If h2o is set to be used instead of the `randomForest` package.
  y <- "Ecosystem" # variable to be predicted
  x <- setdiff(names(meta_comb), c(y, "Site")) # features / predictors
  h2o.init(max_mem_size = h2o_mem_allotted) # adjust as necessary, you might not need so much memory if you don't have a large grid.
  
  eDNA_total.h2o <- as.h2o(meta_comb) # Load data into H2o
  
  if (use_h2o_grid_search == TRUE) { # H2O grid search of parameter space
    eDNA_hypergrid <- list(
      #ntrees = seq(400,2000, by = 800),
      mtries = c(40,100,200,400,600,700),
      sample_rate = c(0.6325, 0.70, 0.80),
      min_rows = seq(1,5, by = 2),
      max_depth = seq(10,40, by = 10)#,
      #nbins = seq(10,20, by = 10)
    )
    
    eDNA_grid_xval <- h2o.grid(
      algorithm = "randomForest",
      grid_id = "eDNArf_grid_xval",
      ntrees = 1000,
      x = x,
      y = y,
      training_frame = eDNA_total.h2o,
      nfolds = 5, # k-fold cross-validation, k = 5
      #fold_column = "Site", # cross-validation by sites, not working or estimated runtime << actual runtime
      hyper_params = eDNA_hypergrid,
      search_criteria = list(strategy = "Cartesian"),
      balance_classes = TRUE # attempts to account for fewer cases in some classes (e.g. dock)
    )
    
    eDNA_grid_performance <- h2o.getGrid( #remember to change the grid ID or make new code for the xval grid 
      grid_id = "eDNArf_grid_xval",
      sort_by = "logloss"
    )
    
    model1 <- h2o.getModel(eDNA_grid_performance@model_ids[[1]]) # Best model, as determined by minimum logloss.
    model1_perf <- h2o.performance(model = model1)
    
  } else { # Using the parameters from the best model from a previous grid search.
    model1 <- h2o.randomForest(
      x = x,
      y = y,
      training_frame = eDNA_total.h2o,
      nfolds = 5,
      ntrees = 1000,
      max_depth = 20,
      min_rows = 3,
      sample_rate = 0.6325,
      mtries = 700,
      keep_cross_validation_predictions = TRUE,
      keep_cross_validation_fold_assignment = TRUE
    )
    model1_perf <- h2o.performance(model = model1)
    model1_path <- h2o.saveModel(object = model1, path = "output/")
    xval_tib <- as_tibble(h2o.cross_validation_holdout_predictions(model1)) # put xval predictions into tibble
    xval_tib$Ecosystem <- meta_comb$Ecosystem # add on actual value for Ecosystem variable
    model1_pred_counts <- with(xval_tib, table(Ecosystem, predict))
    model1_pred_percentage <- model1_pred_counts / rowSums(model1_pred_counts)
    #model1_pred_heatmap <- plot_ly(x = ~ rownames(model1_pred_percentage), y = ~ colnames(model1_pred_percentage), z = ~ model1_pred_percentage, type = "heatmap", colors = rev(lacroix_palettes$Tangerine[1,]))
    model1_pred_heatmap <- plot_ly(x = ~ colnames(model1_pred_percentage), y = ~ rownames(model1_pred_percentage), z = ~ model1_pred_percentage, type = "heatmap", colors = viridis::viridis(255, option = "A"), name = "Frequency", colorbar = list(title = list(text = "Frequency"))) %>%
      layout(
        xaxis = list(
          title = "Predicted"
        ),
        yaxis = list(
          title = "Source"
        )
      )
    model1_pred_heatmap
    ## Figure 7 ---------------------------
    if (use_plotly_cloud == TRUE) {
      api_create(model1_pred_heatmap, filename = "bocas/edna/rf/xval_confusion_heatmap", sharing = plotly_sharing_setting)
    }
  
  var_importances <- h2o.varimp(model1) %>%
    filter(variable != "Site") %>% # temporary. This line to be removed once models are corrected to remove SITE.
    as_tibble %>% # Appending taxonomy information
    mutate(
      Superkingdom = as.character(tax[var_importances$variable,"Superkingdom"]),
      Kingdom = as.character(tax[var_importances$variable,"Kingdom"]), 
      Phylum = as.character(tax[var_importances$variable,"Phylum"]),
      Class = as.character(tax[var_importances$variable,"Class"]),
      Order = as.character(tax[var_importances$variable,"Order"]),
      Family = as.character(tax[var_importances$variable,"Family"]),
      Genus = as.character(tax[var_importances$variable,"Genus"]),
      Species = as.character(tax[var_importances$variable,"species"])
    )
  write_tsv(var_importances, path = metazoan_taxa_importance_file_location)
  }
}


# Technical replicates vs intra-site distances -------------------
# Calculate Bray-Curtis distances between technical replicate samples and compare them to samples from the same site & ecosystem (i.e. subsite, such as all Salt Creek mangroves or all STRI seagrass water samples).

triplicated_samples <- sample_data(curated_bocas_edna.ps) %>% as.data.frame() %>% 
  filter(Sample == "Repeat" ) %>% # sample data of repeat and misc samples. # removed & Ecosystem != "Dock"
  arrange(Site.Code)
triplicated_sites <- triplicated_samples$Site.Code %>% unique()
triplicates_ps <- curated_bocas_edna.ps %>% 
  subset_samples(Site.Code %in% triplicated_sites) %>%
  filter_taxa( function (x) sum(x) > 0, TRUE)

triplicated_sample_data <- sample_data(curated_bocas_edna.ps)
triplicate_distmat <- triplicates_ps %>%
  vegan_otu() %>%
  vegdist() %>%
  as.matrix

triplicated_sample_names <- triplicates_ps %>%
  sample_data() %>%
  as.data.frame() %>%
  rownames_to_column(var = "identifier") %>%
  arrange(Site.Code) %>%
  filter(!is.na(identifier))

triplicate_distmat_sorted <- triplicate_distmat[triplicated_sample_names$identifier, triplicated_sample_names$identifier]

triplicate_data <- tibble()
for (i in triplicated_sites) {
  ids <- triplicated_sample_names[which(triplicated_sample_names$Site.Code == i), "identifier"]
  rep_distmat <- triplicate_distmat_sorted[ids,ids]
  ret_tibble <- tibble(site = i, sampleA = c(base::colnames(rep_distmat)[c(2,3,3)]), sampleB = c(row.names(rep_distmat)[c(1,1,2)]), dist = c(rep_distmat[1,2],rep_distmat[1,3],rep_distmat[2,3]), comparison = "replicate")
  triplicate_data <- rbind(triplicate_data, ret_tibble)
}

# Within & between site variation
sample_data(curated_bocas_edna.ps)$sitehab <- str_remove_all(sample_data(curated_bocas_edna.ps)$site_code3,'[0-9]')
primary_samples_data <- sample_data(curated_bocas_edna.ps) %>% 
  as.data.frame()
primary_samples_data$identifier <- rownames(primary_samples_data)
primary_samples_data <- primary_samples_data %>%
  filter(Fraction == "eDNA" & Sample == "Primary" ) %>% # sample data of repeat and misc samples. #removed & Ecosystem != "Dock"
  arrange(Site.Code)


all_sitehabs <- primary_samples_data$sitehab %>%
  unique()

wi_site_ps <- curated_bocas_edna.ps %>% 
  subset_samples(sitehab %in% all_sitehabs) %>%
  filter_taxa( function (x) sum(x) > 0, TRUE)

wi_site_vegdist <- wi_site_ps %>%
  vegan_otu() %>%
  vegdist() %>%
  as.matrix()

wi_distmat_sorted <- wi_site_vegdist[primary_samples_data$identifier, primary_samples_data$identifier]

intrasite_data <- tibble()
for (i in all_sitehabs) {
  ids <- primary_samples_data[which(primary_samples_data$sitehab == i), "identifier"]
  rep_distmat <- wi_distmat_sorted[ids,ids]
  ret_tibble <- tibble(site = i, sampleA = NA, sampleB = NA, dist = rep_distmat[upper.tri(rep_distmat)], comparison = "intrasite")
  intrasite_data <- rbind(intrasite_data, ret_tibble)
}

trip_comp_tb <- rbind(triplicate_data, intrasite_data) %>%
  mutate(hab_short = str_extract(site,"(?<=-)[CMS]a?"))

hab_map_tb <- tibble(short = c("C", "M", "S", "Sa"), long = c("Coral reef", "Mangrove", "Seagrass", "Unvegetated")) %>%
  column_to_rownames(var = "short")
trip_comp_tb$habitat <- hab_map_tb[trip_comp_tb$hab_short,"long"]
trip_comp_box <- plot_ly(trip_comp_tb, x = ~ habitat, color = ~ comparison, y = ~dist, type = "box") %>%
  layout(boxmode = "group", xaxis = list(title = "Source Ecosystem"), yaxis = list(title = "Bray-Curtis distance", rangemode = "tozero"))

trip_comp_box
if (use_plotly_cloud == TRUE) {
  api_create(trip_comp_box, filename = "bocas/edna/replicates/boxplot", sharing = plotly_sharing_setting)
}


# Ordinations --------------------

## Reef Life Surveys
RLS_nmds <- vegan::metaMDS(fish_table,distance = "bray", k = 3, trymin = 100, trymax = 800, binary = FALSE)
RLS_nmds_2d <- vegan::metaMDS(fish_table,distance = "bray", k = 2, trymin = 100, trymax = 800, binary = FALSE)

RLS_nmds_df <- cbind(RLS_nmds$points,fish_site_metadata_sub)
RLS_nmds_2ddf <- cbind(RLS_nmds_2d$points,fish_site_metadata_sub)

RLS_nmds_p <- plot_ly(type = 'scatter3d', mode = 'markers', data =RLS_nmds_df, x = ~MDS1, y = ~MDS2, z = ~MDS3,
                      text = ~Site, color = ~Ecosystem, symbol = ~BayRegion, 
                      marker = list(opacity = 0.75), colors = "Set1", symbols = c('circle', 'circle-open','x', 'circle-open', 'square')) %>%
  layout(title= paste("RLS Bray-Curtis NMDS, Stress:", round(RLS_nmds$stress, digits = 4)))

RLS_nmds_2d_p <- plot_ly(type = 'scatter', mode = 'markers', data = as.data.frame(RLS_nmds_2ddf), x = ~MDS1, y = ~MDS2,
                         text = ~Site.Name, color = ~Ecosystem, symbol = ~BayRegion, 
                         marker = list(opacity = 0.75, size = 16), colors = "Set1", symbols = c('circle', 'diamond-open-dot','x', 'square')) %>%
  layout(title= paste("RLS Bray-Curtis NMDS, Stress:", round(RLS_nmds_2d$stress, digits = 4)))

RLS_nmds_p
RLS_nmds_2d_p

if (use_plotly_cloud == TRUE) {
  api_create(RLS_nmds_p, filename = "bocas/edna/ordi/RLS_bc_nmds", sharing = plotly_sharing_setting)
  api_create(RLS_nmds_2d_p, filename = "bocas/edna/ordi/RLS_bc_nmds_2d", sharing = plotly_sharing_setting)
}

#RLS PCoA

RLS_pcoa <- cmdscale(vegdist(fish_table, method = "bray", binary = FALSE), k = 2, eig = TRUE)
RLS_pcoa_df <- tibble(Sample = rownames(RLS_pcoa$points), x = RLS_pcoa$points[,1], y = RLS_pcoa$points[,2])
RLS_pcoa_df <- bind_cols(RLS_pcoa_df, fish_site_metadata_sub)
levels(RLS_pcoa_df$Ecosystem)[3] <- "Coral reef"
RLS_pcoa_df$Habitat <- RLS_pcoa_df$Ecosystem


RLS_pcoa_plot <- ggplot(data = RLS_pcoa_df) +
  geom_point(mapping = aes(x = x, y = y, color = Habitat), size = 4, alpha = 0.70) +
  #viridis::scale_color_viridis() + # sorry, but the LaCroix color palettes are too much fun.
  #scale_color_gradientn(colors = lacroix_palette("PeachPear")) +
  #geom_label_repel(mapping = aes(x = x, y = y, label = label, color = year), force = 5) +
  xlab(paste0("PCoA axis 1: ", round(RLS_pcoa$eig[1] / sum(RLS_pcoa$eig) * 100, digits = 1), "%")) +
  ylab(paste0("PCoA axis 2: ", round(RLS_pcoa$eig[2] / sum(RLS_pcoa$eig) * 100, digits = 1), "%")) +
  ggtitle("a") +
  scale_color_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Open Sans" # This is a freely available Google font that I imported into R with the 'extrafont' library.
    ),
    plot.title = element_text(
      face="bold",
      size = 30,
      vjust = -3,
      hjust = 0.05)
  )
RLS_pcoa_plot

ggsave(filename = "output/RLS_pcoa.png", plot = RLS_pcoa_plot, dpi = "retina", width = 5, height = 5, units = "in")



## Fish eDNA
eDNA_fish_vegan <- vegan_otu(eDNA_fish.ps)
eDNA_fish_vegan <- eDNA_fish_vegan[rowSums(eDNA_fish_vegan) > 0,] #remove sites with 0 fish reads
eDNA_fish_bc <- vegdist(eDNA_fish_vegan, method = "bray")
#eDNA_fish_ord <- metaMDS(eDNA_fish_bc, k = 2, trymin = 100, trymax = 800, binary = FALSE)
eDNA_fish_pcoa <- cmdscale(eDNA_fish_bc, k = 2, eig = TRUE)
eDNA_fish_pcoa_df <- tibble(Sample = rownames(eDNA_fish_pcoa$points), x = eDNA_fish_pcoa$points[,1], y = eDNA_fish_pcoa$points[,2])
eDNA_fish_pcoa_df$Ecosystem <- sample_data_renamed$Ecosystem[match(eDNA_fish_pcoa_df$Sample,rownames(sample_data_renamed))]
eDNA_fish_pcoa_df$BayRegion <- sample_data_renamed$BayRegion[match(eDNA_fish_pcoa_df$Sample,rownames(sample_data_renamed))]
eDNA_fish_pcoa_df$Habitat <- eDNA_fish_pcoa_df$Ecosystem

eDNA_fish_pcoa_plot <- ggplot(data = eDNA_fish_pcoa_df) +
  geom_point(mapping = aes(x = x, y = y, color = Habitat), size = 4, alpha = 0.70) +
  #viridis::scale_color_viridis() + # sorry, but the LaCroix color palettes are too much fun.
  #scale_color_gradientn(colors = lacroix_palette("PeachPear")) +
  #geom_label_repel(mapping = aes(x = x, y = y, label = label, color = year), force = 5) +
  xlab(paste0("PCoA axis 1: ", round(eDNA_fish_pcoa$eig[1] / sum(eDNA_fish_pcoa$eig) * 100, digits = 1), "%")) +
  ylab(paste0("PCoA axis 2: ", round(eDNA_fish_pcoa$eig[2] / sum(eDNA_fish_pcoa$eig) * 100, digits = 1), "%")) +
  ggtitle("b") +
  scale_color_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Open Sans" # This is a freely available Google font that I imported into R with the 'extrafont' library.
    ),
    plot.title = element_text(
      face="bold",
      size = 30,
      vjust = -3,
      hjust = 0.05)
  )
eDNA_fish_pcoa_plot
ggsave(eDNA_fish_pcoa_plot, filename = "output/eDNA_fish_PCoA.png", dpi = "retina", width = 5, height = 5, units = "in")

## Metazoan eDNA
curated_edna_metadata_tab <- as_tibble(sample_data(meta_curated_bocas_edna.ps.tr))

meta_OTUs <- vegan_otu(meta_curated_bocas_edna.ps.tr)
meta_bc <- vegdist(meta_OTUs, distance = "bray", binary = FALSE)
meta_pcoa <- cmdscale(meta_bc, k = 2, eig = TRUE)
meta_nmds <- vegan::metaMDS(meta_OTUs, distance = "bray", k = 3, trymin = 100, trymax = 800, binary = FALSE)
meta_nmds_2d <- vegan::metaMDS(meta_OTUs, distance = "bray", k = 2, trymin = 100, trymax = 800, binary = FALSE)
curated_meta_nmds_df <- cbind(meta_nmds$points, curated_edna_metadata_tab)
curated_meta_nmds_2ddf <- cbind(meta_nmds_2d$points, curated_edna_metadata_tab)

meta_nmds_p <- plot_ly(type = 'scatter3d', mode = 'markers', data = curated_meta_nmds_df, x = ~MDS1, y = ~MDS2, z = ~MDS3,
                      text = ~Site, color = ~Ecosystem, symbol = ~BayRegion, 
                      marker = list(opacity = 0.75), colors = "Set1", symbols = c('circle', 'square','x', 'diamond-open-dot')) %>%
  layout(title= paste("Metazoan eDNA Bray-Curtis NMDS, Stress:", round(meta_nmds$stress, digits = 4)))

meta_nmds_p
if (use_plotly_cloud == TRUE) {
  api_create(meta_nmds_p, filename = "bocas/edna/ordi/metazoan_bc_nmds", sharing = plotly_sharing_setting)
}


meta_nmds_2d_p <- plot_ly(type = 'scatter', mode = 'markers', data = curated_meta_nmds_2ddf, x = ~MDS1, y = ~MDS2,
                         text = ~Site, color = ~Ecosystem, symbol = ~BayRegion, 
                         marker = list(opacity = 0.75, size = 16), colors = "Set1", symbols = c('circle', 'square','x', 'diamond-open-dot')) %>%
  layout(title= paste("Metazoan eDNA Bray-Curtis NMDS, Stress:", round(meta_nmds_2d$stress, digits = 4)))

meta_nmds_2d_p
if (use_plotly_cloud == TRUE) {
  api_create(meta_nmds_2d_p, filename = "bocas/edna/ordi/metazoan_bc_nmds_2d", sharing = plotly_sharing_setting)
}


# eDNA metazoan PCoA
meta_pcoa <- cmdscale(meta_bc, k = 2, eig = TRUE)
meta_pcoa_df <- tibble(Sample = rownames(meta_pcoa$points), x = meta_pcoa$points[,1], y = meta_pcoa$points[,2])
meta_pcoa_df <- bind_cols(meta_pcoa_df, curated_edna_metadata_tab)
meta_pcoa_df$Habitat <- meta_pcoa_df$Ecosystem

meta_pcoa_plot <- ggplot(data = meta_pcoa_df) +
  geom_point(mapping = aes(x = x, y = y, color = Habitat), size = 4, alpha = 0.70) +
  #viridis::scale_color_viridis() + # sorry, but the LaCroix color palettes are too much fun.
  #scale_color_gradientn(colors = lacroix_palette("PeachPear")) +
  #geom_label_repel(mapping = aes(x = x, y = y, label = label, color = year), force = 5) +
  xlab(paste0("PCoA axis 1: ", round(meta_pcoa$eig[1] / sum(meta_pcoa$eig) * 100, digits = 1), "%")) +
  ylab(paste0("PCoA axis 2: ", round(meta_pcoa$eig[2] / sum(meta_pcoa$eig) * 100, digits = 1), "%")) +
  ggtitle("c") +
  scale_color_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(
    text = element_text(
      family = "Open Sans" # This is a freely available Google font that I imported into R with the 'extrafont' library.
    ),
    plot.title = element_text(
      face="bold",
      size = 30,
      vjust = -3,
      hjust = 0.05)
  )
meta_pcoa_plot

ggsave(filename = "output/meta_pcoa.png", plot = meta_pcoa_plot, dpi = "retina", width = 5, height = 5, units = "in")


# CCAs #### IN-PROGRESS ######

std_env_var <- read.csv("data/RLS_metadata_updated.csv") %>% #Load standardized metadata file from Elaine
  filter( Site != "PPR") # Site PPR has no sponge AFDW
std_env_var$mangrove_AGB_Z = scale(std_env_var$mangrove_AG_biomass, center = TRUE, scale = TRUE)
  
HAS_sheet <- readxl::read_excel("data/HAS_07072017.xlsx") %>% 
  dplyr::rename(site_code = `Site Abbrev.`, HAS_score = `HAS score`) %>%
  group_by(site_code) %>%
  summarize(mean_HAS = mean(HAS_score))
HAS_sheet$HAS_score_z <- as.vector(scale(HAS_sheet$mean_HAS, center = TRUE, scale = TRUE))
HAS_sheet$HAS_score_z


mangrove_env <- std_env_var %>%
  filter(Habitat == "M") %>%
  column_to_rownames(var = "Subsite")

seagrass_env <- std_env_var %>%
  filter(Habitat == "S") %>%
  column_to_rownames(var = "Subsite")

mangrove_env$HAS_score_z <- HAS_sheet$HAS_score_z[match(mangrove_env$Site, HAS_sheet$site_code)]
mangrove_env_trimmed <- mangrove_env[,c(5,7,8)]

mangrove_bocas_edna <- meta_curated_bocas_edna.ps.tr %>%
  subset_samples(Ecosystem == "Mangrove" & Site %in% mangrove_env$Site) %>%
  merge_samples("Site") %>%
  #subset_samples(site_code3 != "SCR-M") %>% #remove outlier
  filter_taxa( function (x) sum(x) > 0, TRUE) %>%
  transform_sample_counts(function(x) x / sum(x)) %>% 
  filter_taxa( filterfun(kOverA(3,5e-5)), TRUE)

mangrove_otu_for_cca <- vegan_otu(mangrove_bocas_edna)
mangrove_cca <- cca(mangrove_otu_for_cca, mangrove_env_trimmed)
plot(mangrove_cca, scaling = 2, main = "Mangrove eDNA (metazoan) normal CCA", xlab = paste0("CCA1 [", round(mangrove_cca$CCA$eig[1]/sum(mangrove_cca$CCA$eig) * 100, digits = 2), "%]"), ylab = paste0("CCA2 [", round(mangrove_cca$CCA$eig[2]/sum(mangrove_cca$CCA$eig) * 100, digits = 2), "%]"))
vegan:::text.cca(mangrove_cca)


# PERMANOVAs --------------------------------------------

# RLS fish survey

RLS_permanova_nostrata <- adonis2(formula = fish_table ~ BayRegion * Ecosystem + Diver + Site * Ecosystem, data = fish_site_metadata_sub, method = "bray")
RLS_permanova_strata <- adonis2(formula = fish_table ~ BayRegion / Ecosystem, strata = Site, data = fish_site_metadata_sub, method = "bray", strata = "BayRegion")
#adonis2(formula = fish_table ~ Location * Habitat + Diver + Site * Habitat, data = fish_site_metadata_sub, method = "raup", binary = TRUE)

# Plotting code example.
RLS_permanova_nostrata_p <- plotly_permanova(RLS_permanova_nostrata, name = "Fish survey") %>%
  layout(title = "RLS Diver survey (fish)")

RLS_permanova_nostrata_p
if (use_plotly_cloud == TRUE) {
  api_create(RLS_permanova_nostrata_p, filename = "bocas/edna/permanova/RLS_nostrata", sharing = plotly_sharing_setting)
}

# Metazoan eDNA OTUs
# `BayRegion` in the metazoan eDNA data is the same as `Location`
metazoan_permanova_strata <- adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem, strata = Site, data = curated_edna_metadata_tab, method = "bray")
metazoan_permanova_nostrata <- adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem + Site * Ecosystem, data = curated_edna_metadata_tab, method = "bray")
#adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem + Site * Ecosystem, data = curated_edna_metadata_tab, method = "raup", binary = TRUE)

metazoan_permanova_nostrata_p <- plotly_permanova(metazoan_permanova_nostrata, name = "Metazoan eDNA") %>%
  layout(title = "Metazoan eDNA")
metazoan_permanova_nostrata_p

if (use_plotly_cloud == TRUE) {
  api_create(metazoan_permanova_nostrata_p, filename = "bocas/edna/permanova/metazoan_nostrata", sharing = plotly_sharing_setting)
}


permanova_combi_p <- subplot(RLS_permanova_nostrata_p, metazoan_permanova_nostrata_p, shareY = TRUE, shareX = TRUE) %>%
  layout(title = "")
permanova_combi_p

if (use_plotly_cloud == TRUE) {
  api_create(permanova_combi_p, filename = "bocas/edna/permanova/combined_nostrata", sharing = plotly_sharing_setting)
}

# Species list cross-checking ---------------------------


if(file.exists(processed_bocasDB_file)) { # Load processed Bocas spp db
  bocas_metazoans <- read_tsv(processed_bocasDB_file)
} else {
  bocas_metazoans <- read_tsv(unprocessed_bocasDB_file) %>% select(Class, Order, Family, genspec = `Genus/Species`, `Common Name`) # check to make sure quotation marks or other offending chars are removed from the input data or else things will go horribly wrong
  tsns <- get_tsn(bocas_metazoans$genspec)
  accepted_names <- itis_acceptname(tsns[!is.na(tsns)])
  bocas_metazoans$tsn <- tsns
  bocas_metazoans[!is.na(tsns),"accepted_tsn"] <- accepted_names$acceptedtsn
  write_tsv(bocas_metazoans, processed_bocasDB_file)
}


curated_edna.ps.spp <- curated_bocas_edna.ps %>%
  subset_samples(Fraction == "eDNA" & Sample == "Primary") %>% #filter out control, negative and mock samples # if desired, we can also use  `& Ecosystem != "Dock"` to remove docks
  merge_samples("site_code3", fun = mode) %>% # Doesn't matter what the `fun` argument is set to, I think, we overwrite it in the next step
  subset_taxa(Kingdom == "Metazoa") %>% #just used to retrieve a list of species-level IDs
  tax_glom(taxrank = "species") %>%
  filter_taxa(function (x) sum(x) > 0, TRUE)

# How many eDNA OTUs identified down to species level and how many reads does this account for?
sum(sample_sums(curated_edna.ps.spp)) #526898 reads
curated_edna.ps.spp # 387 taxa


edna_species_names <- tax_table(curated_edna.ps.spp)@.Data[,"species"]

if(file.exists(eDNA_spp_file)) {
  edna_species_tibble <- read_tsv(eDNA_spp_file) #Read in file if already existing
} else{
  edna_species_tsns <- get_tsn(edna_species_names)
  edna_species_accepted <- itis_acceptname(edna_species_tsns[!is.na(edna_species_tsns)])
  edna_species_tibble <- tibble(phylum = tax_table(curated_edna.ps.spp)@.Data[,"Phylum"], class = tax_table(curated_edna.ps.spp)@.Data[,"Class"], order = tax_table(curated_edna.ps.spp)@.Data[,"Order"],family = tax_table(curated_edna.ps.spp)@.Data[,"Family"], genus = tax_table(curated_edna.ps.spp)@.Data[,"Genus"], species = edna_species_names, tsn = edna_species_tsns, accepted_tsn = NA)
  edna_species_tibble$in_bocas_db_raw = edna_species_tibble$species %in% bocas_metazoans$genspec
  edna_species_tibble[!is.na(edna_species_tsns),"accepted_tsn"] <- edna_species_accepted$acceptedtsn
  edna_species_tibble$in_bocas_db_acctsn = if_else(!is.na(edna_species_tibble$accepted_tsn), edna_species_tibble$accepted_tsn %in% bocas_metazoans$accepted_tsn, NA)
  edna_species_tibble$acctsn_or_raw_in_bdb <- edna_species_tibble$in_bocas_db_acctsn | edna_species_tibble$in_bocas_db_raw
  write_tsv(edna_species_tibble, eDNA_spp_file)
}



# Load in Global Invasive Species database query
GISD <- read_delim(GISD_data_location, delim = ";")

GISD[which(GISD$Species %in% edna_species_tibble$species),]
edna_species_tibble[which(edna_species_tibble$species %in% GISD$Species),]

summary(curated_edna.ps.spp %>% get_sample(names(which(edna_species_names == "Mus musculus"))) > 0) # Example of how to inspect how many samples have Mus musculus counts



# Checking RLS data against Bocas spp db, GISD ----------
# IN-PROGRESS, INCOMPLETE
fish_check <- fish_tax %>% 
  as_tibble() %>%
  filter(!is.na(species)) %>%
  select(species) %>%
  unique()


fish_check_tsns <- get_tsn(fish_check$species, ask = TRUE) # This should be run interactively when possible so it can take user input on ambiguous species names (e.g. )

fish_check_accepted <- itis_acceptname(fish_check_tsns[!is.na(fish_check_tsns)])

fish_species_tibble <- tibble(species = fish_check$species, tsn = fish_check_tsns, accepted_tsn = NA)
fish_species_tibble$in_bocas_db_raw = fish_species_tibble$species %in% bocas_metazoans$genspec
fish_species_tibble[!is.na(fish_check_tsns),"accepted_tsn"] <- fish_check_accepted$acceptedtsn

fish_species_tibble$in_bocas_db_acctsn = if_else(!is.na(fish_species_tibble$accepted_tsn), fish_species_tibble$accepted_tsn %in% bocas_metazoans$accepted_tsn, NA)
fish_species_tibble$acctsn_or_raw_in_bdb <- fish_species_tibble$in_bocas_db_acctsn | fish_species_tibble$in_bocas_db_raw
# NOTE: Scarus iseri is written as Scarus iserti in the Bocas spp DB which is not an accepted name and is also not on ITIS, so this will need to be manually accounted for at this time.


# Differentially abundant taxa testing ----------

### Fish ###########################
# Figure 4a, RLS fish, mangrove vs seagrass, top 20 diff ---------------
RLS.MS.ps = subset_samples(fish_ps, Ecosystem %in% c("Seagrass","Mangrove")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)
RLS_MvS <- deseq_wrapper(RLS.MS.ps, ~Ecosystem)

RLS_MvS_results <- results(RLS_MvS, contrast = c("Ecosystem","Mangrove","Seagrass"))

fish_lookup <- as.data.frame(fish_code_dictionary) %>% column_to_rownames(var = "Code")
RLS_MvS_results@rownames <- fish_lookup[RLS_MvS_results@rownames, "Species"]
res_to_tsv(RLS_MvS_results, alpha = 1, filename = "output/RLS_MvS.tsv")

# Figure 4b, eDNA fish, mangrove vs seagrass, top 20 diff -------------
eDNA_fish.ps <- curated_bocas_edna.ps %>%
  subset_samples(Fraction == "eDNA" & Sample == "Primary") %>% #filter out control, negative and mock samples # if desired, we can also use  `& Ecosystem != "Dock"` to remove docks
  merge_samples("site_code3", fun = mode) %>% # Doesn't matter what the `fun` argument is set to, I think, we overwrite it in the next step
  subset_taxa(Class == "Actinopteri") %>% # Actinopteri only
  filter_taxa(function (x) sum(x) > 0, TRUE)# Remove taxa that are no longer represented because of subsetting

sample_data(eDNA_fish.ps) <- sample_data(sample_data_renamed)

eDNA_fish.MS.ps <- eDNA_fish.ps %>%
  subset_samples(Ecosystem %in% c("Seagrass","Mangrove")) %>% # Metazoan only
  filter_taxa(function (x) sum(x) > 0, TRUE)

eDNA_fish_MvS <- deseq_wrapper(eDNA_fish.MS.ps, formula = ~ Ecosystem)

eDNA_fish_MvS.var <- getVarianceStabilizedData(eDNA_fish_MvS, blind = FALSE) # Variance-stabilized data. Not yet used
#eDNA_fish_MvS = DESeq(eDNA_fish_MvS, fitType="local")

eDNA_fish_MvS_results <- results(eDNA_fish_MvS, contrast = c("Ecosystem","Mangrove","Seagrass"))
#eDNA_fish_MvS_results@rownames <- eDNA_fish_MvS_results@rownames.... # Replace SHA1 names with lowest taxonomic name
res_to_tsv(eDNA_fish_MvS_results, alpha = 1, filename = "output/eDNA_fish_MvS.tsv") # alpha = 1 means show all results regardless of p-adj



# Figure "S1 4a", RLS fish, reef vs sand, top 20 diff ---------------------
RLS.RU.ps = subset_samples(fish_ps, Ecosystem %in% c("Reef","Sand")) # Subset to your contrast or else the size scaling factor calculation will incorporate all samples (w/unnecessary habitats)

deseq_wrapper(RLS.RU.ps, formula = ~ Ecosystem)


RLS_RvU_results <- results(RLS_RvU, contrast = c("Ecosystem","Coral reef","Unvegetated"))

fish_lookup <- as.data.frame(fish_code_dictionary) %>% column_to_rownames(var = "Code")
RLS_RvU_results@rownames <- fish_lookup[RLS_RvU_results@rownames, "Species"]
res_to_tsv(RLS_RvU_results, alpha = 1, filename = "output/RLS_RvU.tsv")



# Figure "S1 4b", eDNA fish, reef vs sand, top 20 diff ------------------

eDNA_fish.RU.ps <- eDNA_fish.ps %>%
  subset_samples(Ecosystem %in% c("Coral reef","Unvegetated")) %>%
  filter_taxa(function (x) sum(x) > 0, TRUE)


eDNA_fish_RvU <- deseq_wrapper(eDNA_fish.RU.ps, ~Ecosystem)
eDNA_fish_RvU.var <- getVarianceStabilizedData(eDNA_fish_RvU, blind = FALSE) # Variance-stabilized data. Not yet used
#eDNA_fish_RvU = DESeq(eDNA_fish_RvU, fitType="local")

eDNA_fish_RvU_results <- results(eDNA_fish_RvU, contrast = c("Ecosystem","Coral reef","Unvegetated"))
#eDNA_fish_RvU_results@rownames <- eDNA_fish_RvU_results@rownames.... # Replace SHA1 names with lowest taxonomic name
res_to_tsv(eDNA_fish_RvU_results, alpha = 1, filename = "output/eDNA_fish_RvU.tsv") # alpha = 1 means show all results regardless of p-adj


# Figure 4A/b set 2 inner vs outer bay, RLS (a) and eDNA (b)  -------------
RLS.bay.ps = subset_samples(fish_ps, BayRegion %in% c("Inner Bay","Outer Bay")) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

RLS_bay <- deseq_wrapper(RLS.bay.ps, formula = ~BayRegion)

#Size factor and dispersion estimates
RLS_bay.var <- getVarianceStabilizedData(RLS_bay, blind = FALSE) # Variance-stabilized data. Not yet used
#RLS_bay = DESeq(RLS_bay, fitType="local")

RLS_bay <- DESeq(RLS_bay, test = "Wald", fitType = "parametric")

RLS_bay_results <- results(RLS_bay, contrast = c("BayRegion","Outer Bay","Inner Bay"))
#RLS_bay_results@rownames <- RLS_bay_results@rownames.... # Replace SHA1 names with lowest taxonomic name
res_to_tsv(RLS_bay_results, alpha = 1, filename = "output/RLS_bay.tsv") # alpha = 1 means show all results regardless of p-adj

eDNA_fish.bay.ps <- eDNA_fish.ps %>%
  subset_samples(BayRegion %in% c("Inner Bay","Outer Bay")) %>%
  filter_taxa(function (x) sum(x) > 0, TRUE)


eDNA_fish_bay <- deseq_wrapper(eDNA_fish.bay.ps, formula = ~BayRegion)
eDNA_fish_bay.var <- getVarianceStabilizedData(eDNA_fish_bay, blind = FALSE) # Variance-stabilized data. Not yet used
#eDNA_fish_bay = DESeq(eDNA_fish_bay, fitType="local")

eDNA_fish_bay <- DESeq(eDNA_fish_bay, test = "Wald", fitType = "parametric")

eDNA_fish_bay_results <- results(eDNA_fish_bay, contrast = c("BayRegion","Outer Bay","Inner Bay"))
#eDNA_fish_bay_results@rownames <- eDNA_fish_bay_results@rownames.... # Replace SHA1 names with lowest taxonomic name
res_to_tsv(eDNA_fish_bay_results, alpha = 1, filename = "output/eDNA_fish_bay.tsv") # alpha = 1 means show all results regardless of p-adj


### Metazoa ########################

# Figure 5a, eDNA metazoa, mangrove vs seagrass


eDNA_meta.MS.ps <- eDNA_meta.ps %>%
  subset_samples(Ecosystem %in% c("Seagrass","Mangrove")) %>% # Metazoan only
  filter_taxa(function (x) sum(x) > 0, TRUE)

eDNA_meta_MvS <- deseq_wrapper(eDNA_meta.MS.ps, formula = ~Ecosystem)
eDNA_meta_MvS.var <- varianceStabilizingTransformation(eDNA_meta_MvS, blind = FALSE) # Variance-stabilized version
#plotUpDownSigTaxa(eDNA_meta_MvS_results, rld = eDNA_meta_MvS.var, metadata_cols = "Ecosystem", title = "metazoan eDNA, mangrove vs seagrass")

eDNA_meta_MvS_results <- results(eDNA_meta_MvS, contrast = c("Ecosystem","Mangrove","Seagrass"))
#eDNA_fish_MvS_results@rownames <- eDNA_fish_MvS_results@rownames.... # Replace SHA1 names with lowest taxonomic name
res_to_tsv(eDNA_meta_MvS_results, alpha = 1, filename = "output/eDNA_meta_MvS.tsv") # alpha = 1 means show all results regardless of p-adj

plotUpDownSigTaxa_eco(results = eDNA_meta_MvS_results, rld = eDNA_meta_MvS.var, metadata_cols = "Ecosystem", ntaxa = 40, title = "eDNA metazoa sorted")

# Figure 5b, eDNA metazoa, inner vs outer bay -------------------
 eDNA_meta.bay.ps <-eDNA_meta.ps %>%
  subset_samples(BayRegion %in% c("Inner Bay", "Outer Bay")) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE)

eDNA_meta_bay <- deseq_wrapper(eDNA_meta.bay.ps, formula = ~ BayRegion)
eDNA_meta_bay.var <- getVarianceStabilizedData(eDNA_meta_bay, blind = FALSE)
eDNA_meta_bay_results <- results(eDNA_meta_bay)
res_to_tsv(eDNA_meta_bay_results, alpha = 1, filename = "output/eDNA_meta_bay.tsv")
