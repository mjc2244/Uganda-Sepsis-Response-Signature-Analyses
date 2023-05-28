#Clear R environment
rm(list = ls())

#Import the raw gene counts
raw <- read.csv(file.choose(), header=TRUE, check.names = FALSE)

#Convert all gene names fixed row identifier
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

#Stratify patients using 19 gene SRS
library(SepstratifieR)

predictions19 <- stratifyPatients(loggeneset19, gene_set = "extended")

plotAlignedSamples(predictions19)

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

#Stratify patients using 7 gene SRS
predictions7 <- stratifyPatients(loggeneset7, gene_set = "davenport")

plotAlignedSamples(predictions7)

#Merge SRS7 predictions/assignments to clinical datasets 
#Quantitiative SRS
srsq7 <- predictions7@SRSq
srsq7

#Qualitative
srs7 <- predictions7@SRS
srs7

#Merge SRS assignments 
reservesrs3 <- merge(reservesrs2, srs7, by='row.names',all=TRUE)
rownames(reservesrs3) = reservesrs3$Row.names
reservesrs3$Row.names = NULL
colnames(reservesrs3)[291] = "srs7"

reservesrs4 <- merge(reservesrs3, srsq7, by='row.names',all=TRUE)
rownames(reservesrs4) = reservesrs4$Row.names
reservesrs4$Row.names = NULL
colnames(reservesrs4)[292] = "srsq7"

#Ok now lets make the plots showing relationship between SRSq, MEWS, UVA
library(ggplot2)
library(ggprism)
library(mgcv)
library(mgcViz)

#MEWS - 19 gene
mews_19 <- gam(MEWS_score ~ s(srsq19), data = reservesrs4, method = "REML")
summary(mews_19)

mews_19 <- getViz(mews_19)
check(mews_19)

mews_19_plot <- plot(mews_19, trans = function(x){ 
  mews_19$family$linkinv(coef(mews_19)[1] + x)
}, seWithMean = TRUE, rug = TRUE, shade = TRUE)
mews_19_plot <- mews_19_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
mews_19_plot <- mews_19_plot + scale_x_continuous(name ="19-gene SRSq") + scale_y_continuous(name = "MEWS") + theme_bw()
mews_19_plot <- mews_19_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
mews_19_plot <- mews_19_plot + 
  annotate("text", x = 0.4, y = 5, label = "p < 0.001",hjust=0.5, vjust=0, 
           size= 5, col="black") 
mews_19_plot

#UVA - 19 gene
uva_19 <- gam(UVA_score ~ s(srsq19), data = reservesrs4, method = "REML")
summary(uva_19)

uva_19 <- getViz(uva_19)
check(uva_19)

uva_19_plot <- plot(uva_19, trans = function(x){ 
  uva_19$family$linkinv(coef(uva_19)[1] + x)
}, seWithMean = TRUE, rug = TRUE, shade = TRUE)
uva_19_plot <- uva_19_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
uva_19_plot <- uva_19_plot + scale_x_continuous(name ="19-gene SRSq") + scale_y_continuous(name = "UVA score") + theme_bw()
uva_19_plot <- uva_19_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
uva_19_plot <- uva_19_plot + 
  annotate("text", x = 0.4, y = 5, label = "p = 0.003",hjust=0.5, vjust=0, 
           size= 5, col="black") 
uva_19_plot

#MEWS - 7 gene
mews_7 <- gam(MEWS_score ~ s(srsq7), data = reservesrs4, method = "REML")
summary(mews_7)

mews_7 <- getViz(mews_7)
check(mews_7)

mews_7_plot <- plot(mews_7, trans = function(x){ 
  mews_7$family$linkinv(coef(mews_7)[1] + x)
}, seWithMean = TRUE, rug = TRUE, shade = TRUE)
mews_7_plot <- mews_7_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
mews_7_plot <- mews_7_plot + scale_x_continuous(name ="7-gene SRSq") + scale_y_continuous(name = "MEWS") + theme_bw()
mews_7_plot <- mews_7_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
mews_7_plot <- mews_7_plot + 
  annotate("text", x = 0.25, y = 6, label = "p < 0.001",hjust=0.5, vjust=0, 
           size= 5, col="black") 
mews_7_plot

#UVA - 7 gene
uva_7 <- gam(UVA_score ~ s(srsq7), data = reservesrs4, method = "REML")
summary(uva_7)

uva_7 <- getViz(uva_7)
check(uva_7)

uva_7_plot <- plot(uva_7, trans = function(x){ 
  uva_7$family$linkinv(coef(uva_7)[1] + x)
}, seWithMean = TRUE, rug = TRUE, shade = TRUE)
uva_7_plot <- uva_7_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
uva_7_plot <- uva_7_plot + scale_x_continuous(name ="7-gene SRSq") + scale_y_continuous(name = "UVA score") + theme_bw()
uva_7_plot <- uva_7_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
uva_7_plot <-  uva_7_plot + 
  annotate("text", x = 0.25, y = 5, label = "p = 0.001",hjust=0.5, vjust=0, 
           size= 5, col="black") 
uva_7_plot

#Hospital outcome
hospital_19 <- gam(hospdeathtransf ~ s(srsq19), data = reservesrs4, family=binomial, method = "REML")
summary(hospital_19)

hospital_19 <- getViz(hospital_19)
check(hospital_19)

hospital_19_plot <- plot(hospital_19, trans=plogis, seWithMean = TRUE, rug = TRUE, shade = TRUE)
hospital_19_plot <- hospital_19_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
hospital_19_plot <- hospital_19_plot + scale_x_continuous(name ="19-gene SRSq") + scale_y_continuous(name = "Hospital death or transfer") + theme_bw()
hospital_19_plot <- hospital_19_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
hospital_19_plot <- hospital_19_plot + 
  annotate("text", x = 0.4, y = 0.75, label = "p = 0.439",hjust=0.5, vjust=0, 
           size= 5, col="black") 
hospital_19_plot

#Hospital outcome
hospital_7 <- gam(hospdeathtransf ~ s(srsq7), data = reservesrs4, family=binomial, method = "REML")
summary(hospital_7)

hospital_7 <- getViz(hospital_7)
check(hospital_7)

hospital_7_plot <- plot(hospital_7, trans=plogis, seWithMean = TRUE, rug = TRUE, shade = TRUE)
hospital_7_plot <- hospital_7_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
hospital_7_plot <- hospital_7_plot + scale_x_continuous(name ="7-gene SRSq") + scale_y_continuous(name = "Hospital death or transfer") + theme_bw()
hospital_7_plot <- hospital_7_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
hospital_7_plot <- hospital_7_plot + 
  annotate("text", x = 0.2, y = 0.75, label = "p = 0.378",hjust=0.5, vjust=0, 
           size= 5, col="black") 
hospital_7_plot

#30d mortality
death30d_19 <- gam(death30d ~ s(srsq19), data = reservesrs4, family=binomial, method = "REML")
summary(death30d_19)

death30d_19 <- getViz(death30d_19)
check(death30d_19)

death30d_19_plot <- plot(death30d_19, trans=plogis, seWithMean = TRUE, rug = TRUE, shade = TRUE)
death30d_19_plot <- death30d_19_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
death30d_19_plot <- death30d_19_plot + scale_x_continuous(name ="19-gene SRSq") + scale_y_continuous(name = "Death at 30 days") + theme_bw()
death30d_19_plot <- death30d_19_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
death30d_19_plot <- death30d_19_plot + 
  annotate("text", x = 0.4, y = 0.8, label = "p = 0.002",hjust=0.5, vjust=0, 
           size= 5, col="black") 
death30d_19_plot

#30d mortality
death30d_7 <- gam(death30d ~ s(srsq7), data = reservesrs4, family=binomial, method = "REML")
summary(death30d_7)

death30d_7 <- getViz(death30d_7)
check(death30d_7)

death30d_7_plot <- plot(death30d_7, trans=plogis, seWithMean = TRUE, rug = TRUE, shade = TRUE)
death30d_7_plot <- death30d_7_plot + l_fitLine(colour = "blue", size = 1) + l_rug(mapping = aes(x=x), alpha = 0.8) +
  l_ciLine(level = 0.95, mul = NULL, colour = "gray66", linetype = 2, size =1)
death30d_7_plot <- death30d_7_plot + scale_x_continuous(name ="7-gene SRSq") + scale_y_continuous(name = "Death at 30 days") + theme_bw()
death30d_7_plot <- death30d_7_plot + theme_prism(base_size = 14) + theme(legend.position = "none")
death30d_7_plot <- death30d_7_plot + 
  annotate("text", x = 0.25, y = 0.8, label = "p = 0.001",hjust=0.5, vjust=0, 
           size= 5, col="black") 
death30d_7_plot




