# This is a lulu script for post-clustering of Zotu table created by unoise3

# Enable extracting command-line arguments
args = commandArgs(trailingOnly=TRUE)
minMatch_lulu <- args[1]

require(lulu)
otutab <- read.csv("zotuTable.txt",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table("match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

# Curation step with lulu
curated_result <- lulu(otutab, matchlist, minimum_match = minMatch_lulu) # This runs the default parameter of lulu (i.e. minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95)
# there are more parameters to play with in lulu command. Check lulu paper to understand how they will affect your results 

curated_result$curated_table # Curated OTU table
curated_result$curated_count # Number of OTUs retained
curated_result$curated_otus # IDs of curated OTUs
curated_result$discarded_count # OTUs discarded
curated_result$otu_map # total - total read count, spread - the number of samples the OTU is present in
                       # parent_id - ID of OTU with which this OTU was merged (or self)
                       # curated - ("parent" or "merged"), was this OTU kept as a valid OTU (parent) or merged with another
                       # rank - The rank of the OTU in terms of decreasing spread and read count

curated_result$original_table # Original OTU table

write.table(curated_result$curated_table, "curated_zotuTable.tab", sep="\t")  # write curated result 
write.table(curated_result$otu_map, "lulu_zotu_map.tab", sep="\t")            # write the map info
