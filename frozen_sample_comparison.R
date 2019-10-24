# This script compares and tests DNA extraction and PCR yields for frozen and freshly-filtered water samples.
library(tidyverse)
library(lubridate)
library(plotly) # This requires a plotly account + additional setup for Plotly in R.
BCS10_metadata <- readxl::read_excel("data/BCS_10_COI (eDNA).xlsx", sheet = 1)
BCS12_metadata <- readxl::read_excel("data/BCS_12_COI (Coiba eDNA & repeats).xlsx", sheet = 2) %>%
  select(MLID = `DNA ID`, PCR_yield = `mlCOI/jgHCO [PCR] cleaned undiluted`)

base_metadata <- readr::read_tsv("data/eDNA_sampledata.tsv") %>% # Presumably extracted later means filtered later and not actually DNA extracted.
  filter(!str_detect(`Full Site Name`, "[Nn]eg")) %>% # Remove negatives
  mutate(sample_date = mdy(Sample.Date), extraction_date = mdy(`Date filtered (m/d/y)`)) %>%
  select(MLID, full_name = `Full Site Name`, extraction_yield = `QuBit reading (ng/uL)`, sample_date, extraction_date) %>%
  mutate(extracted_later = extraction_date > sample_date)
repeated_MLIDs <- BCS12_metadata$`DNA ID`[which(BCS12_metadata$`DNA ID` %in% BCS10_metadata$`Matt Code`)]

BCS10_trim <- BCS10_metadata %>%
  select(MLID = `Matt Code`, PCR_yield = `mlCOI/jgHCO [PCR] cleaned undiluted`) %>%
  filter(MLID != "Mock", MLID %in% base_metadata$MLID)

combined <- as_tibble(cbind(BCS10_trim, base_metadata[match(BCS10_trim$MLID, base_metadata$MLID),-1]))
sum(combined$MLID %in% repeated_MLIDs)
length(unique(repeated_MLIDs)) #repeats were not part of the triplicates, so we can just blindly replace them.

combined$PCR_yield[combined$MLID %in% repeated_MLIDs] <- BCS12_metadata$PCR_yield[match(combined$MLID[combined$MLID %in% repeated_MLIDs], BCS12_metadata$MLID)] # Replace yields in our dataframe with the final libraries for repeated libraries.

combined <- combined %>%
  mutate(PCR_yield = as.numeric(PCR_yield), extraction_yield = as.numeric(extraction_yield))

combined2 <- tibble(MLID = combined$MLID, PCR_yield = combined$PCR_yield, extracted_later = combined$extracted_later)

PCR_yield_box <- ggplot(data = combined) +
  geom_boxplot(mapping = aes( y = PCR_yield, x = extracted_later) )
PCR_yield_box

PCR_yield_box <- plot_ly(data = combined, y = ~PCR_yield, x = ~ extracted_later, boxpoints = "all", type = "box") %>%
  layout(
    xaxis = list(
      title = "Frozen sample"
    ),
    yaxis = list(
      title = "PCR yield (ng / µl)"
    )
  )

PCR_yield_box
t.test(combined$PCR_yield[which(combined$extracted_later)], combined$PCR_yield[which(!combined$extracted_later)])

ext_yield_box <- plot_ly(data = combined, y = ~extraction_yield, x = ~ extracted_later, boxpoints = "all", type = "box") %>%
  layout(
    xaxis = list(
      title = "Frozen sample"
    ),
    yaxis = list(
      title = "DNA extraction yield (ng / µl)"
    )
  )

ext_yield_box
t.test(combined$extraction_yield[which(combined$extracted_later)], combined$extraction_yield[which(!combined$extracted_later)])

api_create(PCR_yield_box, filename = "bocas/edna/frozen_samples/PCR_yield", sharing = "secret")
api_create(ext_yield_box, filename = "bocas/edna/frozen_samples/extraction_yield", sharing = "secret")
