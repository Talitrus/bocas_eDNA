#!/usr/bin/env Rscript

# pass an integer (number of cores) to the Rscript to use multiple cores. Otherwise it defaults to 1 core.
# Used only for running on the Smithsonian HPC Hydra or for other clusters where the number of cores cannot be reliably estimated using `nproc`.
# Currently, there are no parallel operations in this script, but I have this in many of scripts for convenience in case I need it later.
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  nslots <- 1
  print("No argument provided. Defaulting to single core processing.")
} else if (length(args)==1 & is.integer(as.integer(args[1]))) {
  # default output file
  nslots <- as.integer(args[1])
  print(paste("Using", nslots, "cores."))
} else {
  stop("More than one argument supplied. You should only supply one argument, which should be the number of cores allotted as an integer.")
}

# Load libraries ------------------------

library(phyloseq)
library(vegan)
library(tidyverse)
library(parallel)
library(plotly)
library(lubridate)
library(randomForest)
library(genefilter) # for the kOverA function
library(taxize)
library(rsample)

# Flags & settings variables ----------------------

use_plotly_cloud <- FALSE # set this to TRUE to enable the use of `api_create` to make/update Plot.ly cloud figures.
plotly_sharing_setting <- "secret" # Plot.ly cloud privacy setting for figures generated, not used if `use_plotly_cloud` is FALSE.

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

# Functions --------------
vegan_otu <- function(physeq) { #convert phyloseq OTU table into vegan OTU matrix
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

filter_taxa_to_other <- function(physeq, filterFunc, merged_taxon_name = "Other taxa") { #returns phyloseq object
  taxa_to_merge <- (!filter_taxa(physeq, filterFunc, FALSE)) # named logical vector where samples NOT meeting the filter threshold are TRUE
  taxa_to_merge_names <- which(taxa_to_merge == TRUE) %>% names() #get the names
  new_physeq <- merge_taxa(physeq, taxa_to_merge_names, archetype = 1)
  tax_table(new_physeq)[taxa_to_merge_names[1],] <- replicate(ncol(tax_table(new_physeq)), merged_taxon_name) # merges with first taxon, so we can replace the lineage information for the first taxon to whatever we want displayed
  return(new_physeq)
}



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

phyglom_curated_bocas_edna.ps.tr <- curated_bocas_edna.ps %>%
  subset_samples(Fraction == "eDNA") %>% #filter out control, negative and mock samples, NOTE: currently includes non-primary samples
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  filter_taxa_to_other(function(x) max(x) > 0.03) # All taxa that are never > 0.03 rel. abun. merged into "Other taxa" 

# Load fish data -------------------------------------

fish_data <- read_csv(RLS_survey_location) %>%
  filter(Diver != 'Daniella Heflin') %>% # We had observations from a third diver for some sites—removed for consistency
  select(Line_ID, Diver, Buddy, Site.No., event, Site.Name, Habitat, Latitude, Longitude, Date, vis, Direction, Time, Depth, Method, Block, Code, Species, Common.name, `Total Abundance`, `Total Biomass`) %>%
  mutate(Code = str_to_lower(Code)) %>% # Code is an abbrev. for species, lowercased for consistency
  unique() # One row is an exact duplicate, likely an error. Removing duplicates
fish_data$Species[fish_data$Species == "Clupeid spp."] <- "Clupeidae" # manual correction

if(file.exists(fish_site_metadata_location)) { #If a metadata file exists, read it in
  # Currently, this is used because I manually added a column called Location to the metadata file, which corresponds to "BayRegion" in the eDNA
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
    select(event, Site.Name,Habitat, Latitude, Longitude, Depth, Method, Block, datetime, Diver) %>%
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


fish_ps <- phyloseq(otu_table(as.matrix(fish_table), taxa_are_rows = FALSE), tax_table(fish_tax), sample_data(fish_site_metadata_sub))

# Taxonomic composition -------------------------------------------
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

curated_edna_meta.tr.f <- meta_curated_bocas_edna.ps.tr %>%
  filter_taxa(filterfun(kOverA(5,5e-5)), TRUE) #filter for taxa with approximately at least 5 reads per sample in at least 5 samples)
selected_env <- sample_data(curated_edna_meta.tr.f) %>%
  select(Ecosystem)
meta_comb <- cbind(otu_table(curated_edna_meta.tr.f),selected_env) #combined dataframe with OTUs plus the Ecosystem column.
meta_rf <- randomForest(meta_comb[,1:(ncol(meta_comb)-1)], meta_comb[,ncol(meta_comb)], ntree = 2000, importance = TRUE, proximity = TRUE, mtry = 200)
raw_importance <- importance(meta_rf) %>% as.data.frame() %>% rownames_to_column(var = "SHA1") %>% as_tibble() # Feature importance

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
write_tsv(meta_imp_table, path = metazoan_taxa_importance_file_location)

# Trying out training + testing data split for cross validation.
eDNA_split <- initial_split(meta_comb, prop = 0.7, strata = 'Ecosystem')
eDNA_training <- training(eDNA_split)
eDNA_testing <- testing(eDNA_split)

meta_rf_cv <- randomForest(x = eDNA_training[,1:(ncol(eDNA_training)-1)], y = eDNA_training[,ncol(eDNA_training)], xtest = eDNA_testing[,1:(ncol(eDNA_testing)-1)], ytest = eDNA_testing[,ncol(eDNA_testing)], ntree = 2000, importance = TRUE, proximity = TRUE, mtry = 200)

# The test set error is relatively similar to the OOB error.


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
  layout(boxmode = "group", xaxis = list(title = "Habitat"), yaxis = list(title = "Bray-Curtis distance", rangemode = "tozero"))

trip_comp_box
if (use_plotly_cloud == TRUE) {
  api_create(trip_comp_box, filename = "bocas/edna/replicates/boxplot", sharing = plotly_sharing_setting)
}

#Need to make a site_code4 column or something that's like CCR-Sa with no number, then loop through site/habitat combos in the distance matrix using the primary_samples_data and stash those distances into a data frame.

# Ordinations --------------------

## Reef Life Surveys
RLS_nmds <- vegan::metaMDS(fish_table,distance = "bray", k = 3, trymin = 100, trymax = 800, binary = FALSE)
RLS_nmds_2d <- vegan::metaMDS(fish_table,distance = "bray", k = 2, trymin = 100, trymax = 800, binary = FALSE)

RLS_nmds_df <- cbind(RLS_nmds$points,fish_site_metadata_sub)
RLS_nmds_2ddf <- cbind(RLS_nmds_2d$points,fish_site_metadata_sub)

RLS_nmds_p <- plot_ly(type = 'scatter3d', mode = 'markers', data =RLS_nmds_df, x = ~MDS1, y = ~MDS2, z = ~MDS3,
                      text = ~Site, color = ~Habitat, symbol = ~Location, 
                      marker = list(opacity = 0.75), colors = "Set1", symbols = c('circle', 'circle-open','x', 'circle-open', 'square')) %>%
  layout(title= paste("RLS Bray-Curtis NMDS, Stress:", round(RLS_nmds$stress, digits = 4)))

RLS_nmds_2d_p <- plot_ly(type = 'scatter', mode = 'markers', data = as.data.frame(RLS_nmds_2ddf), x = ~MDS1, y = ~MDS2,
                         text = ~Site.Name, color = ~Habitat, symbol = ~Location, 
                         marker = list(opacity = 0.75, size = 16), colors = "Set1", symbols = c('circle', 'diamond-open-dot','x', 'square')) %>%
  layout(title= paste("RLS Bray-Curtis NMDS, Stress:", round(RLS_nmds_2d$stress, digits = 4)))

RLS_nmds_p

if (use_plotly_cloud == TRUE) {
  api_create(RLS_nmds_p, filename = "bocas/edna/ordi/RLS_bc_nmds", sharing = plotly_sharing_setting)
  api_create(RLS_nmds_2d_p, filename = "bocas/edna/ordi/RLS_bc_nmds_2d", sharing = plotly_sharing_setting)
}



## Metazoan eDNA
curated_edna_metadata_tab <- as_tibble(sample_data(meta_curated_bocas_edna.ps.tr))

meta_OTUs <- vegan_otu(meta_curated_bocas_edna.ps.tr)
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

if (use_plotly_cloud == TRUE) {
  api_create(meta_nmds_2d_p, filename = "bocas/edna/ordi/metazoan_bc_nmds_2d", sharing = plotly_sharing_setting)
}





# PERMANOVAs --------------------------------------------

# RLS fish survey

RLS_permanova_nostrata <- adonis2(formula = fish_table ~ Location * Habitat + Diver + Site * Habitat, data = fish_site_metadata_sub, method = "bray")
RLS_permanova_strata <- adonis2(formula = fish_table ~ Location / Habitat, strata = Site, data = fish_site_metadata_sub, method = "bray", strata = "Location")

#adonis2(formula = fish_table ~ Location * Habitat + Diver + Site * Habitat, data = fish_site_metadata_sub, method = "raup", binary = TRUE)

# Metazoan eDNA OTUs
# `BayRegion` in the metazoan eDNA data is the same as `Location`
metazoan_permanova_strata <- adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem, strata = Site, data = curated_edna_metadata_tab, method = "bray")
metazoan_permanova_nostrata <- adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem + Site * Ecosystem, data = curated_edna_metadata_tab, method = "bray")
#adonis2(formula = meta_OTUs ~ BayRegion * Ecosystem + Site * Ecosystem, data = curated_edna_metadata_tab, method = "raup", binary = TRUE)

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


curated_edna.ps.spp <- meta_curated_bocas_edna.ps.tr %>% #just used to retrieve a list of species-level IDs
  tax_glom(taxrank = "species") %>%
  filter_taxa(function (x) sum(x) > 0, TRUE)
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
