#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
# Modified: Jessica Pearce
# 04/11/2022
# Modified: Adam Bennett
# 22/02/2023
#===========================================================================

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(seqinr))

# Define options for command line
option_list = list(
  make_option(c("-l", "--lca_file"), action="store", default=NA, type='character',
              help="Table produced by LCA script (e.g., '09_taxonomy_assigned/lca_taxAssigned_results_qCov100_id97_diff1.tab')"),
  make_option(c("-z", "--zotu_file"), action="store", default=NA, type='character',
              help="OTU or ZOTU table (e.g., '08_lulu/curated_zotuTable.tab')"),
  make_option(c("-f", "--fasta_file"), action="store", default=NA, type='character',
              help="fasta file (e.g., '06_Unique_ZOTUs/*.fasta_zotus.fasta')"),
  make_option(c("-m", "--metadata_file"), action="store", default=NA, type='character',
              help="metadata file in csv format (e.g., 'test_metadata.csv')"),          
  make_option(c("-c", "--cores"), action="store", default=NA, type='integer',
              help="number of cores (e.g., 100)"),
  make_option(c("-p", "--phyloseq_file"), action="store", default=NA, type='character',
              help="The name of the output file"),
  make_option(c("-o", "--optimise"), action="store", default=FALSE, type='logical',
              help="TRUE or FALSE; setting this to TRUE will cause the script take longer to run"))


#........................................................................
# Setup
#........................................................................

opt = parse_args(OptionParser(option_list=option_list))

lca_file      <- opt$lca_file
zotu_file     <- opt$zotu_file
fasta_file    <- opt$fasta_file
metadata_file <- opt$metadata_file
cores         <- opt$cores
phyloseq_file <- opt$phyloseq_file
optimise      <- opt$optimise

lca_tab <- read_tsv(lca_file)
otu_tab   <- read_tsv(zotu_file)

colnames(otu_tab)[1]   <- "ASV"
colnames(lca_tab)[8]   <- "ASV"
fasta <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
  
otu_tab["ASV_sequence"] <- NA
  
for(row in 1:nrow(otu_tab)) {
  ASV <- otu_tab[row, "ASV"]
  otu_tab[row, "ASV_sequence"] <- fasta[[ASV$ASV]][1]
}
  
asv_seq <- otu_tab %>%
  select(ASV, ASV_sequence)

merged_tab <- merge(lca_tab, asv_seq, by = "ASV", all.x = TRUE)

taxa     <- merged_tab
otu      <- merged_tab
seq_tab  <- merged_tab

if (grep(".csv", metadata_file)){
  meta <- read_csv(metadata_file)
} else if (grep(".tsv", metadata_file)) {
  meta <- read_tsv(metadata_file)
}


#........................................................................
# Prepare metadata
#........................................................................

meta             <- as.data.frame(meta)
rownames(meta)   <- meta$sample_id
meta$sample_id <- NULL


#........................................................................
# Prepare taxa data
#........................................................................

taxa["LCA"] <- ""
taxa %>%
  select(ASV, domain, phylum, class, order, family, genus, species, LCA, ASV_sequence) -> taxa
taxa        <- as.data.frame(taxa)

# We need to fill the LCA column starting with 'species', then 'genus', etc
# until we reach a taxonomy level that isn't 'na' or 'dropped'
levels <- c("species", "genus", "family", "order", "class", "phylum", "domain")
for (row in 1:nrow(taxa)) {
  LCA      <- NA
  level_ID <- 1
  while(LCA == "dropped" | is.na(LCA)) {
    LCA      <- taxa[row, levels[level_ID]]
    level_ID <- level_ID + 1
  }
  taxa[row, "LCA"] <- LCA
}

rownames(taxa) <- taxa$ASV
taxa$ASV       <- NULL
taxa           <- as.matrix(taxa)


#........................................................................
# Prepare otu data
#........................................................................

otu           <- as.data.frame(otu)
rownames(otu) <- otu$ASV
otu[,c('ASV', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'Contam', 'ASV_sequence')] <- list(NULL)


#........................................................................
# Create phylogenetic tree
#........................................................................

# Phylogenetic tree code based on code from
# https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
DNA_set = DNAStringSet(seq_tab$ASV_sequence)
names(DNA_set) = paste0(seq_tab$ASV)

alignment = AlignSeqs(DNA_set, anchor=NA, processors=cores)

phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
dm          <- dist.ml(phang_align)
treeNJ      <- NJ(dm)

fit    <- pml(treeNJ, data=phang_align)
fitGTR <- update(fit, k=4, inv=0.2)

if (optimise) {
  # Use 'optim.pml()' to optimise the model parameters
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
}


#........................................................................
# Create phyloseq object
#........................................................................

OTU    <- otu_table(otu, taxa_are_rows = TRUE)
TAX    <- tax_table(taxa)
META   <- sample_data(meta)
TREE   <- phy_tree(fitGTR$tree)

physeq <- phyloseq(OTU, TAX, META, TREE)

saveRDS(physeq, file = paste0(phyloseq_file))