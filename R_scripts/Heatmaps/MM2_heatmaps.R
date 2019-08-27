##
##    MM2_heatmaps.R
##
##
##    DESCRIPTION OF OPERATION FLOW:
##      Part 1: Load required packages
##      Part 2: Source the MM2 samples and DE model file
##      Part 3: PCA Plot
##      Part 4: Heatmaps
##


## PART 1: Load required packages

# Load required packages for this file
library(gplots)

## PART 2: Load .txt file with mapped read count and remove genes + samples with no valuable data
source("R_scripts/shared/Datasets_MM2.R")

## PART 3: PCA Plot
# For comparative and descriptive analysis we use the rlog (see manual for descriptions). 
# blind=FALSE to keep sample group info for calculating normalizations. 
rldMM2 <- rlog(dds, blind = FALSE)

# Perform PCA of samples
pdf(file="output/pdf/PCA/190827_PCA_MM2_all.pdf", height = 24, width = 6)
plotPCA(rldMM2, intgroup=c("condition"))
dev.off()


## PART 4A: Heatmap all samples 
# Extract the rld transformed values
rldMM2_values <- assay(rldMM2)


## Work on ordered heatmaps here
order_to_sort = c("PG22.P20", "PG23.P20", "PG24.P20", "PG22.P21", "PG24.P21", "2PG23.P21",
                  "PG22.P13", "PG24.P13", "PG26.P13", "PG24.P14", "PG26.P14", "2PG22.P14",
                  "PG22.P10", "PG23.P10", "PG24.P10", "PG26.P10","PG22.P11", "PG23.P11", "PG24.P11", "PG26.P11",
                  "PG22.P12", "PG23.P12", "PG24.P12", "PG26.P12", "PG22.P16", "PG23.P16", "PG24.P16", "PG26.P16")

rldMM2_values = rldMM2_values[, order_to_sort]


# Calculate standard deviations per gene (all 9 samples) and order based on std dev (highest at top)
geneSDdds <- rep(NA, nrow(rldMM2_values))
for (i in 1:nrow(rldMM2_values)) {geneSDdds[i] <- sd(rldMM2_values[i,])}
rldMM2_values_sd_ord <- rldMM2_values[order(-geneSDdds),]

my_palette <- colorRampPalette(c("blue","white", "red"))(n = 30)

col_breaks  = c(seq(4,16, length=31))

my_palette_z_score <- colorRampPalette(c("blue","white", "red"))(n = 25)
col_breaks_z_score = c(seq(-2,2, length=26))

pdf(file="output/pdf/Heatmaps/190827_clustering-MM2_all_Top250.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rldMM2_values_sd_ord[1:250,]),
           col=my_palette_z_score,
           breaks=col_breaks_z_score,
           scale=c("row"),
           trace="none", 
           density.info="none", 
           margins =c(4,4),     
           cexCol=0.8,
           cexRow=0.5,
           keysize=7,
           key=T,
           Rowv=T, 
           Colv=F,
           colsep=c(3,6,9,12,16,20,24),
           dendrogram="none",
           lmat        = rbind(c(0,3),c(2,1),c(0,4)),
           lwid        = c(1.5,1), 
           lhei        = c(1.5,3,0.48) 
)
dev.off()

rm(rldMM2_values_sd_ord)

# PART 4B: GFPneg samples cytokine profiles

# Extract EV-GFPneg samples only 
rldMM2_values <- rldMM2_values[,c(which(grepl("P20",colnames(rldMM2_values))),
                                  which(grepl("P21",colnames(rldMM2_values))),
                                  which(grepl("P13",colnames(rldMM2_values))),
                                  which(grepl("P10",colnames(rldMM2_values))),
                                  which(grepl("P12",colnames(rldMM2_values))))]

## PART 5: Extract data for specific gene groups from the rlog values data obtained in part 4

my_palette  = colorRampPalette(c("blue4","white", "red"))(n = 24)
col_breaks  = c(seq(-2,2, length=25)) 
colsep      = c(0,3,6,9,13,17)
sepwidth    = c(0.05,0.05)
sepcolor    = "white"
cexCol      = 0.6
cexRow      = 0.8
lmat        = rbind(c(0,3),c(2,1),c(0,4))
lwid        = c(1.5,2) # second position for column width
lhei        = c(1.5,2,0.3) # second position for row height
offsetRow   = -15.5

# load groups of genes as manually currated in the Groups_for_analysis.R file
source("R_scripts/shared/Groups_for_analysis.R")


# Now obtain the right DE orders. 
dds_resTP12TP20 <- results(dds, c("condition","TP12","TP20"))
dds_resTP12TP20Ord_log2fold <- dds_resTP12TP20[order(-dds_resTP12TP20$log2FoldChange),]

# Work on the heatmaps for the different cytokine heatmaps

# Xue_IFNg_mouse graphs

dds_resTP12TP20Ord_log2fold_IFNg = subset(dds_resTP12TP20Ord_log2fold, rownames(dds_resTP12TP20Ord_log2fold) %in% Xue_IFNg_mouse)
IFNg_order_TP12TP20 = rownames(dds_resTP12TP20Ord_log2fold_IFNg)

rldMM2_values_IFNg_resTP12TP20 = subset(rldMM2_values, rownames(rldMM2_values) %in% Xue_IFNg_mouse)
rldMM2_values_IFNg_resTP12TP20 = as.data.frame(rldMM2_values_IFNg_resTP12TP20)

rldMM2_values_IFNg_resTP12TP20 = rldMM2_values_IFNg_resTP12TP20[match(IFNg_order_TP12TP20, rownames(rldMM2_values_IFNg_resTP12TP20)),]

pdf(file="output/pdf/Heatmaps/190827_mouse_IFNg_GFP_neg_TP12TP20.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rldMM2_values_IFNg_resTP12TP20),
           main="mouse_IFNg", 
           col=my_palette,
           breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rldMM2_values_IFNg_resTP12TP20)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", 
           density.info="none", 
           margins =c(25,12),     
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=offsetRow,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()

# Xue_IL4_mouse graphs

# Set up the orders for TP12TP20 comparison
dds_resTP12TP20Ord_log2fold_IL4 = subset(dds_resTP12TP20Ord_log2fold, rownames(dds_resTP12TP20Ord_log2fold) %in% Xue_IL4_mouse)
IL4_order_TP12TP20 = rownames(dds_resTP12TP20Ord_log2fold_IL4)

rldMM2_values_IL4_resTP12TP20 = subset(rldMM2_values, rownames(rldMM2_values) %in% Xue_IL4_mouse)
rldMM2_values_IL4_resTP12TP20 = as.data.frame(rldMM2_values_IL4_resTP12TP20)

rldMM2_values_IL4_resTP12TP20 = rldMM2_values_IL4_resTP12TP20[match(IL4_order_TP12TP20, rownames(rldMM2_values_IL4_resTP12TP20)),]

pdf(file="output/pdf/Heatmaps/190827_mouse_IL4_GFP_neg_TP12TP20.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rldMM2_values_IL4_resTP12TP20),
           main="mouse_IL4",
           col=my_palette,
           breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rldMM2_values_IL4_resTP12TP20)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", 
           density.info="none", 
           margins =c(25,12),     
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=offsetRow,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()


# Xue_IL10_mouse graphs

# Set up the orders for TP12TP20 comparison
dds_resTP12TP20Ord_log2fold_IL10 = subset(dds_resTP12TP20Ord_log2fold, rownames(dds_resTP12TP20Ord_log2fold) %in% Xue_IL10_mouse)
IL10_order_TP12TP20 = rownames(dds_resTP12TP20Ord_log2fold_IL10)

rldMM2_values_IL10_resTP12TP20 = subset(rldMM2_values, rownames(rldMM2_values) %in% Xue_IL10_mouse)
rldMM2_values_IL10_resTP12TP20 = as.data.frame(rldMM2_values_IL10_resTP12TP20)

rldMM2_values_IL10_resTP12TP20 = rldMM2_values_IL10_resTP12TP20[match(IL10_order_TP12TP20, rownames(rldMM2_values_IL10_resTP12TP20)),]

rm(dds_resTP12TP20Ord_log2fold_IL10, IL10_order_TP12TP20)

pdf(file="output/pdf/Heatmaps/190827_mouse_IL10_GFP_neg_TP12TP20.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rldMM2_values_IL10_resTP12TP20),
           main="mouse_IL10", 
           col=my_palette,
           breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rldMM2_values_IL10_resTP12TP20)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", 
           density.info="none",
           margins =c(25,12),     
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=offsetRow,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()


# GSEA_mouse_IL6_STAT3 graphs

# Set up the orders for TP12TP20 comparison
dds_resTP12TP20Ord_log2fold_IL6_STAT3 = subset(dds_resTP12TP20Ord_log2fold, rownames(dds_resTP12TP20Ord_log2fold) %in% GSEA_mouse_IL6_STAT3)
IL6_STAT3_order_TP12TP20 = rownames(dds_resTP12TP20Ord_log2fold_IL6_STAT3)

rldMM2_values_IL6_STAT3_resTP12TP20 = subset(rldMM2_values, rownames(rldMM2_values) %in% GSEA_mouse_IL6_STAT3)
rldMM2_values_IL6_STAT3_resTP12TP20 = as.data.frame(rldMM2_values_IL6_STAT3_resTP12TP20)

rldMM2_values_IL6_STAT3_resTP12TP20 = rldMM2_values_IL6_STAT3_resTP12TP20[match(IL6_STAT3_order_TP12TP20, rownames(rldMM2_values_IL6_STAT3_resTP12TP20)),]

rm(dds_resTP12TP20Ord_log2fold_IL6_STAT3, IL6_STAT3_order_TP12TP20)

pdf(file="output/pdf/Heatmaps/190827_mouse_IL6_GFP_neg_TP12TP20.pdf", height = 24, width = 6)
heatmap.2 (as.matrix(rldMM2_values_IL6_STAT3_resTP12TP20),
           main="mouse_IL6_STAT3", 
           col=my_palette,
           breaks=col_breaks,
           colsep=colsep,
           rowsep=c(1:nrow(rldMM2_values_IL6_STAT3_resTP12TP20)),
           sepwidth=sepwidth,
           sepcolor=sepcolor, 
           trace="none", 
           density.info="none", 
           margins =c(25,12),     
           cexCol=cexCol,
           cexRow=cexRow,
           scale=c("row"),
           key=T,
           Colv=F,
           Rowv=F,
           dendrogram="none",
           adjRow=c(1,0.5),
           offsetRow=offsetRow,
           lmat = lmat,
           lwid = lwid,
           lhei = lhei
)
dev.off()