#Clear R environment
rm(list = ls())

#Import the raw gene counts
raw <- read.csv(file.choose(), header=TRUE, check.names = FALSE)

#Convert all gene names to fixed row identifier
rownames(raw) = raw$pid
raw$pid = NULL

#CPM transformation
library(scTenifoldNet)
rawcpm <- cpmNormalization(raw)

#Select genes required to eval 19-gene SRS signature
library(dplyr)
data_keep_rows19 <- c("SLC25A38", 
                    "DNAJA3", 
                    "NAT10", 
                    "THOC1", 
                    "MRPS9",
                    "PGS1", 
                    "UBAP1",
                    "USP5",
                    "TTC3",
                    "SH3GLB1",
                    "BMS1",
                    "FBXO31",
                    "ARL14EP",
                    "CCNB1IP1",
                    "DYRK2",
                    "ADGRE3",
                    "MDC1",
                    "TDRD9",
                    "ZAP70") 
geneset19 <- rawcpm[rownames(rawcpm) %in% data_keep_rows19, ]

#Select genes required to eval 7-gene SRS signature
data_keep_rows7 <- c("ARL14EP", 
                     "CCNB1IP1", 
                     "DYRK2", 
                     "ADGRE3", 
                     "MDC1", 
                     "TDRD9", 
                     "ZAP70") 
geneset7 <- rawcpm[rownames(rawcpm) %in% data_keep_rows7, ]

#Transpose datasets
geneset19 <- as.matrix(geneset19)
geneset19 <- t(geneset19)
geneset19 <- as.data.frame(geneset19)

geneset7 <- as.matrix(geneset7)
geneset7 <- t(geneset7)
geneset7 <- as.data.frame(geneset7)

#Log2-transform CPM-normalized gene counts
loggeneset19 <- log2(geneset19)
loggeneset7 <- log2(geneset7)

#Rename column names to Ensembl as required by SepstratifieR
loggeneset19 <- loggeneset19 %>% 
  rename(
    ENSG00000144659 =	SLC25A38,
    ENSG00000103423 =	DNAJA3,
    ENSG00000135372 =	NAT10,
    ENSG00000079134 =	THOC1,
    ENSG00000135972 =	MRPS9,
    ENSG00000087157 =	PGS1,
    ENSG00000165006 =	UBAP1,
    ENSG00000111667 =	USP5,
    ENSG00000182670 =	TTC3,
    ENSG00000097033 =	SH3GLB1,
    ENSG00000165733 =	BMS1,
    ENSG00000103264 =	FBXO31,
    ENSG00000152219 =	ARL14EP,
    ENSG00000100814 =	CCNB1IP1,
    ENSG00000127334 =	DYRK2,
    ENSG00000131355 =	ADGRE3,
    ENSG00000137337 =	MDC1,
    ENSG00000156414 =	TDRD9,
    ENSG00000115085 =	ZAP70)

loggeneset7 <- loggeneset7 %>% 
  rename(
    ENSG00000152219 = ARL14EP,
    ENSG00000100814 = CCNB1IP1,
    ENSG00000127334 =	DYRK2,
    ENSG00000131355 =	ADGRE3,
    ENSG00000137337 =	MDC1,
    ENSG00000156414 =	TDRD9,
    ENSG00000115085 =	ZAP70)

#Stratify patients using 19 gene SRS, check for outliers, sensitivity analysis of k 
library(SepstratifieR)

predictions19 <- stratifyPatients(loggeneset19, gene_set = "extended")

plotAlignedSamples(predictions19)

sensitivity_results19 <- runSensitivityAnalysis(loggeneset19, gene_set = "extended")

plotAlignedSamples(predictions19, pcs=c(1,2), color_by = "mNN_outlier")
plotAlignedSamples(predictions19, pcs=c(1,3), color_by = "mNN_outlier")

#Merge SRS19 predictions/assignments to clinical datasets 
#Quantitiative SRS
srsq19 <- predictions19@SRSq
srsq19

#Qualitative
srs19 <- predictions19@SRS
srs19

#Import clinical data, row_id to fixed, for those with RNAseq data
reserve <- read.csv(file.choose(), header=TRUE)
rownames(reserve) = reserve$pid
reserve$pid = NULL

#Merge w/ SRS assignments 
reservesrs <- merge(reserve, srs19, by='row.names',all=TRUE)
rownames(reservesrs) = reservesrs$Row.names
reservesrs$Row.names = NULL
colnames(reservesrs)[289] = "srs19"

reservesrs2 <- merge(reservesrs, srsq19, by='row.names',all=TRUE)
rownames(reservesrs2) = reservesrs2$Row.names
reservesrs2$Row.names = NULL
colnames(reservesrs2)[290] = "srsq19"

#Stratify patients using 7 gene SRS, check for outliers, sensitivity analysis of k 
predictions7 <- stratifyPatients(loggeneset7, gene_set = "davenport")

plotAlignedSamples(predictions7)

sensitivity_results7 <- runSensitivityAnalysis(loggeneset7, gene_set = "davenport")

plotAlignedSamples(predictions7, pcs=c(1,2), color_by = "mNN_outlier")
plotAlignedSamples(predictions7, pcs=c(1,3), color_by = "mNN_outlier")

#Merge SRS7 predictions/assignments to clinical datasets 
#Quantitiative SRS
srsq7 <- predictions7@SRSq
srsq7

#Qualitative
srs7 <- predictions7@SRS
srs7

#Merge SRS assignments. 
reservesrs3 <- merge(reservesrs2, srs7, by='row.names',all=TRUE)
rownames(reservesrs3) = reservesrs3$Row.names
reservesrs3$Row.names = NULL
colnames(reservesrs3)[291] = "srs7"

reservesrs4 <- merge(reservesrs3, srsq7, by='row.names',all=TRUE)
rownames(reservesrs4) = reservesrs4$Row.names
reservesrs4$Row.names = NULL
colnames(reservesrs4)[292] = "srsq7"

