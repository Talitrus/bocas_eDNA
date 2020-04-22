# Environmental DNA survey captures patterns of fish and invertebrate diversity across a tropical seascape
Published article:
Nguyen, B.N., Shen, E.W., Seemann, J. et al. Environmental DNA survey captures patterns of fish and invertebrate diversity across a tropical seascape. Sci Rep 10, 6729 (2020). https://doi.org/10.1038/s41598-020-63565-9

## Authors
Bryan N Nguyen†, Elaine W Shen†, Janina Seemann, Adrienne MS Correa, James L O’Donnell, Andrew H Altieri, Nancy Knowlton, Keith A Crandall, Scott P Egan, W Owen McMillan, Matthieu Leray

† These authors contributed equally.

## Introduction
This GitHub page holds the code for processing the 
Bocas del Toro 2017 environmental DNA project starting after trimming and demultiplexing with Flexbar, amplicon sequence variant-calling with [DADA2](https://benjjneb.github.io/dada2/), 
clustering with [VSEARCH](https://github.com/torognes/vsearch), post-OTU curation with [LULU](https://github.com/tobiasgf/lulu), and taxonomic assignment with [BLCA](https://github.com/qunfengdong/BLCA) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Because these steps above were run on a high-performance cluster and the scripts and organization of other clusters and bioinformaticians may vary, we provide the parameters used for the above softwares below instead of the exact scripts used.

### Flexbar

Flexbar (version 3.0.3) was run twice to remove primer indices from both the forward and reverse reads. Reads had to have identical primer indices on both ends to pass quality control. Flexbar was run with the following parameters:

```bash
flexbar -r <read 1 FASTQ> -p <read 2 FASTQ> --adapter-trim-end ANY -a <TruSeq adapters FASTA> -ao 7 -b <primer index FASTA> --barcode-trim-end LTAIL
```
Any files under 100k in size were removed for quality control then Flexbar was run again on the output files from the first Flexbar run with the following parameters:

```bash
flexbar -r <read 2 FASTQ> -p <read 1 FASTQ> -b <primer index FASTA> --barcode-trim-end LTAIL
```

### DADA2

After adapter-removal and demultiplexing with Flexbar, reads were processed with DADA2 (version 1.9.0) using the following parameters:

```R
BCS_out <- filterAndTrim(forward_reads, forward_reads_filtered, reverse_reads, reverse_reads_filtered, rm.phix = TRUE, maxN = 0, maxEE = c(2,2), truncQ = 10, trimLeft = 26, minLen = 100, multithread = parallel::detectCores())
```

Sequence variants were re-oriented into the forward direction by reverse-complementing if the reverse complement had a better Hamming distance than the original sequence after Needleman-Wunsch alignment to a reference sequence from a 313-bp section of the COI gene from a _Tedania ignis_ sequence, GenBank accession number DQ133904.1. This is to correct for the tendency of the adapter-ligation library preparation method to cause random read directions.

```
>DQ133904.1 subset
ATTATCTGGGATTCAGGCTCATTCCGGGGGTTCGGTAGATTTGGTTATTTTTAGTTTACATTTAGCGGGTATTTCTTCTATATTGGCGGCTATGAATTTTATAACCACTATTATTAATATGAGGGCACCAGGGATAACAATGGATAGAACGCCATTGTTTGTTTGGTCAATTTTAGTAACTGCGGTTTTATTATTATTATCTTTACCAGTATTAGCAGGCGCAATTACTATGTTATTAACGGATAGAAATTTTAATACTGCTTTTTTTGATCCAGCAGGTGGAGGAGACCCGATTTTATATCAACATTTATTT
```

Chimeras were then detected and removed using DADA2's `removeBimeraDenovo()` function.

### VSEARCH OTU clustering
Sequence variants from DADA2 were checked for chimeras again with VSEARCH and clustered at 97% identity.

### OTU curation with LULU
OTUs from VSEARCH were then curated using LULU with the following parameters: `minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95`.
### Taxonomic assignment

#### BLCA
Sequence variants from DADA2 were assigned taxonomy using a [slightly modified version of BLCA](https://github.com/Talitrus/BLCA) forked from commit `fb2bd12`. Sequences were compared against the MIDORI UNIQUE database version 20180221. Hits needed to have a minimum identity score of 70 and a minimum coverage of 0.75 to be used.

#### Iterative BLAST search
OTUs that remained unidentified with BLCA were compared to the whole NCBI NT database (retrieved May 2018) using BLAST searches (word size = 7; max e-value = 5e-13) and assigned the taxonomic information of the lowest common ancestor of the top 100 hits.

After this, OTU and taxonomy tables were combined into a single phyloseq object in R (`data/curated_Bocas_eDNA_phyloseq.rds` in this GitHub repository).

## Organization
This GitHub repository contains the code for the R analyses accompanying this project, downstream from the analyses described above. The file `edna_manuscript_minimal.R` contains the bulk of the primary analyses. The file `differential_abundance_heatmaps.R` contains the code for differential abundance analyses and heatmap figure generation. The file `frozen_sample_comparison.R` contains code for testing for differences in DNA extraction yield and PCR yield between frozen and freshly-filtered water samples.

These scripts should run properly if you clone the repo and set the repo as the working directory in R when you run the R scripts.
