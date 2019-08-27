##
##    MA_plot_data_MM2_samples.R
##
##

## PART 1: Load required packages

# Load required packages
library(gplots)

## PART 2: source the file generating the dataframe. We use the minus get microglia datasets file
source("R_scripts/shared/Datasets_MM2.R")

# TP14 vs. TP13
dds_resTP14TP13 <- results(dds, c("condition","TP14","TP13"))
DESeq2::plotMA(dds_resTP14TP13,  alpha = 0.05, main="MA plot TP14 vs TP13", ylim=c(-10,10))

# TP11 vs. TP10
dds_resTP11TP10 <- results(dds, c("condition","TP11","TP10"))
DESeq2::plotMA(dds_resTP11TP10,  alpha = 0.05, main="MA plot TP11 vs TP10", ylim=c(-10,10))

# TP16 vs. TP12
dds_resTP16TP12 <- results(dds, c("condition","TP16","TP12"))
DESeq2::plotMA(dds_resTP16TP12,  alpha = 0.05, main="MA plot TP16 vs TP12", ylim=c(-10,10))

# TP20 vs. TP21
dds_resTP20TP21 <- results(dds, c("condition","TP20","TP21"))
DESeq2::plotMA(dds_resTP20TP21,  alpha = 0.05, main="MA plot TP20 vs TP21", ylim=c(-10,10))

# TP13 vs. TP20
dds_resTP13TP20 <- results(dds, c("condition","TP13","TP20"))
DESeq2::plotMA(dds_resTP13TP20,  alpha = 0.05, main="MA plot TP13 vs TP20", ylim=c(-10,10))
