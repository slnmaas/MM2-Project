library(DESeq2)
library(limma)


# PART 2: Load samples and do basic quality control

cntMM2 <- read.table(file="input/star_genes_erc.counts_NS20_EA.txt")
cntMM2_batch2 <- read.table(file="input/NS26_star_genes_erc.counts-MM2-2.txt")

# Rename samples from batch2 so that we can identify them subsequently
names(cntMM2_batch2)[names(cntMM2_batch2) == 'PG22.P14_S54.gene_erc.countsHT.txt'] <- '2PG22.P14_S54.gene_erc.countsHT.txt'
names(cntMM2_batch2)[names(cntMM2_batch2) == 'PG23.P21_S55.gene_erc.countsHT.txt'] <- '2PG23.P21_S55.gene_erc.countsHT.txt'

# Merge the two sets together
cntAll <- merge(cntMM2, cntMM2_batch2, by="row.names", all=TRUE)

# Set the names of the rows properly again after merging
rownames(cntAll) <- cntAll[,c("Row.names")]
cntAll[,c("Row.names")] <- NULL

# Remove sets we no longer need
rm(cntMM2, cntMM2_batch2)

# Remove features (certain rownames) we no longer need
colnames(cntAll) <- strsplit2(colnames(cntAll),"_")[,1]
htseq_drop_rows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique","ERCC-")
drop <- rowSums(sapply(htseq_drop_rows, grepl, rownames(cntAll)))>0
cntMM2All <- cntAll[!drop,]

# Remove stuff we no longer need
rm(cntAll, drop, htseq_drop_rows)

# Remove colls we don't want to analyze
cntMM2All$BR1 <- NULL
cntMM2All$pmac1 <- NULL
cntMM2All$SHIMG1.3 <- NULL

# Also remove PG22.P14 since we have a better sample with more reads: 2PG22.P14
cntMM2All$PG22.P14 <- NULL

# Print dimensions of 
dim(cntMM2All)

# Keep genes with more than 1 read in at least 2 samples
keep <- rowSums(cntMM2All>=1)>=2
table(keep)
cntMM2All <- cntMM2All[keep,]
dim(cntMM2All)

# Drop low read count (less than 5 reads) samples. First show these
hist(colSums(cntMM2All>5), breaks=20)

# Print colSums output to spot the samples with less counts
colSums(cntMM2All>5)

# Drop samples where less than 6000 genes were detected with more than 5 reads
drop <- colSums(cntMM2All>5)<6000
table(drop)
cntMM2All <- cntMM2All[,!drop]

hist(colSums(cntMM2All>5), breaks=20)

rm(drop, keep)


# PART 2: DESeq2 model
colDataDMM <- data.frame(condition=rep(NA,length(colnames(cntMM2All))), row.names=colnames(cntMM2All))
colDataDMM$condition[grep("P....P10",rownames(colDataDMM))] <- "TP10"
colDataDMM$condition[grep("P....P11",rownames(colDataDMM))] <- "TP11"
colDataDMM$condition[grep("P....P12",rownames(colDataDMM))] <- "TP12"
colDataDMM$condition[grep("P....P16",rownames(colDataDMM))] <- "TP16"
colDataDMM$condition[grep("P....P13",rownames(colDataDMM))] <- "TP13"
colDataDMM$condition[grep("P....P14",rownames(colDataDMM))] <- "TP14"
colDataDMM$condition[grep("P....P20",rownames(colDataDMM))] <- "TP20"
colDataDMM$condition[grep("P....P21",rownames(colDataDMM))] <- "TP21"
dds <- DESeqDataSetFromMatrix(countData = cntMM2All,colData = colDataDMM,design = ~condition)
dds$condition <- factor(dds$condition,levels=c("TP13",
                                               "TP14",
                                               "TP12",
                                               "TP16",
                                               "TP10",
                                               "TP11",
                                               "TP20",
                                               "TP21"))
dds <- DESeq(dds)
rm(colDataDMM, cntMM2All)