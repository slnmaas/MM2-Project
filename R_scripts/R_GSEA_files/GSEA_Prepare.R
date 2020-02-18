##
##    GSEA_Prepare.R
##

## PART 1: source the file generating the dataframe
source("R_scripts/shared/Datasets_MM2.R")

## PART 2: generate the .rnk file based on the DE order TP13 versus TP20
dds_resTP13TP20 = results(dds, c("condition","TP13","TP20"))
dds_resTP13TP20 = dds_resTP13TP20[order(dds_resTP13TP20$log2FoldChange, decreasing = TRUE),]
dds_resTP13TP20 = dds_resTP13TP20[c("log2FoldChange")]
write.table(dds_resTP13TP20,file="output/tables/GSEA/dds_resTP13TP20_Ranked.rnk",col.names=F,quote=F,sep="\t")
