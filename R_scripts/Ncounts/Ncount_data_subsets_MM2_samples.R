##
##    Ncount_data_subsets_MM2_samples.R
##    
##    PART 1: Load required packages
##    PART 2: Source the MM2 samples and DE model file
##    PART 3: Obtain normalized reads for the subset of genes and write to disk
##

## PART 1: Load required packages

# Load required packages
library(gplots)

## PART 2: Source the MM2 samples and DE model file
source("R_scripts/shared/Datasets_MM2.R")

## PART 3: Obtain normalized reads for the subset of genes and write to disk
normalized_counts <- counts(dds, normalized=TRUE)

## Work on ordered set
order_to_sort = c("PG22.P20", "PG23.P20", "PG24.P20", "PG22.P21", "PG24.P21", "2PG23.P21",
                  "PG22.P13", "PG24.P13", "PG26.P13", "PG24.P14", "PG26.P14", "2PG22.P14",
                  "PG22.P10", "PG23.P10", "PG24.P10", "PG26.P10","PG22.P11", "PG23.P11", "PG24.P11", "PG26.P11",
                  "PG22.P12", "PG23.P12", "PG24.P12", "PG26.P12", "PG22.P16", "PG23.P16", "PG24.P16", "PG26.P16")

normalized_counts = normalized_counts[, order_to_sort]

# Vector of genes of interest
Genes_of_interest = c("Csf1r", "Ly6c1", "Ly6c2", "Ccr2", "Cx3cr1", "Spn", "Sell", "Treml4", "Fcgr1", "Mertk", "Itgax",
                      "Ptprc", "Adgre1", "Mrc1", "Ccr7", "Nr4a1", "Cd74", "H2-Aa", "H2-Eb1", "Il1b", "Arg1", "H2-DMb1",
                      "Il4ra", "Tnf", "Mrc1", "Cd80", "Cd86", "Nos2", "Adgre1")
Genes_of_interest = unique(Genes_of_interest)

# Extract ncounts and save file
Ncount_genes_of_interest = subset(normalized_counts, rownames(normalized_counts) %in% Genes_of_interest)
write.table(Ncount_genes_of_interest,file="output/tables/190827_Ncount_genes_of_interest.txt",col.names=T,quote=F,sep="\t")
