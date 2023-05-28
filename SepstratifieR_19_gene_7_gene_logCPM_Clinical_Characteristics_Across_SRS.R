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

#Table of clinical/micro characteristics across SRS groups
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessduration",
              "historyoffever", 
              "nightsweats", 
              "headache", 
              "cough", 
              "diarrhea", 
              "shortnessofbreath",
              "painwithurination", 
              "abxmalarialprior", 
              "temp38plus", 
              "templess36", 
              "heartrate3", 
              "resprate3", 
              "sbp3",
              "o2sat3", 
              "ams", 
              "qsofa2p", 
              "qsofa1p",
              "sirs",
              "MEWS_score",
              "UVA_score", 
              "shock", 
              "srd", 
              "sevanemia",
              "hivrdtresult", 
              "hivstage34", 
              "newhivdx", 
              "artprior",
              "tmpsmxprior", 
              "malariardtresult", 
              "microtbdx", 
              "urinetblamresult",
              "influenzapcrresult", 
              "hospdeathtransf", 
              "daysspentinhospital",
              "kps70less", 
              "death30d",
              "srsq19",
              "trclust")

#Define categorical variables
catVars <- c("sex", 
             "historyoffever", 
             "nightsweats", 
             "headache", 
             "cough", 
             "diarrhea", 
             "shortnessofbreath",
             "painwithurination", 
             "abxmalarialprior", 
             "temp38plus", 
             "templess36", 
             "ams", 
             "qsofa2p", 
             "qsofa1p", 
             "sirs",
             "shock", 
             "srd", 
             "sevanemia",
             "hivrdtresult", 
             "hivstage34", 
             "newhivdx", 
             "artprior",
             "tmpsmxprior", 
             "malariardtresult", 
             "microtbdx", 
             "urinetblamresult",
             "influenzapcrresult", 
             "hospdeathtransf", 
             "kps70less", 
             "death30d",
             "trclust")

tablesrs19 <- CreateTableOne(vars = listVars, 
                             data = reservesrs2, 
                             factorVars = catVars,
                             strata = "srs19",
                             includeNA = TRUE,
                             addOverall = TRUE)
tablesrs19
summary(tablesrs19)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tablesrs19 <- print(tablesrs19, nonnormal = c("age",
                                              "illnessduration",
                                              "heartrate3", 
                                              "resprate3", 
                                              "sbp3",
                                              "o2sat3",
                                              "qsofa_score",
                                              "sirs",
                                              "MEWS_score",
                                              "UVA_score",
                                              "daysspentinhospital",
                                              "srsq19"),
                    exact = c("sex", 
                              "historyoffever", 
                              "nightsweats", 
                              "headache", 
                              "cough", 
                              "diarrhea", 
                              "shortnessofbreath",
                              "painwithurination", 
                              "abxmalarialprior", 
                              "temp38plus", 
                              "templess36", 
                              "ams", 
                              "qsofa2p", 
                              "qsofa1p", 
                              "sirs",
                              "shock", 
                              "srd", 
                              "sevanemia",
                              "hivrdtresult", 
                              "hivstage34", 
                              "newhivdx", 
                              "artprior",
                              "tmpsmxprior", 
                              "malariardtresult", 
                              "microtbdx", 
                              "urinetblamresult",
                              "influenzapcrresult", 
                              "hospdeathtransf", 
                              "kps70less", 
                              "death30d",
                              "trclust"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tablesrs19), file = 'srs19groups.html')

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

#Table of clinical/micro characteristics across SRS groups
library(tableone)
#Create a variable list which we want in Table
listVars <- c("sex", 
              "age",
              "illnessduration",
              "historyoffever", 
              "nightsweats", 
              "headache", 
              "cough", 
              "diarrhea", 
              "shortnessofbreath",
              "painwithurination", 
              "abxmalarialprior", 
              "temp38plus", 
              "templess36", 
              "heartrate3", 
              "resprate3", 
              "sbp3",
              "o2sat3", 
              "ams", 
              "qsofa2p", 
              "qsofa1p",
              "sirs",
              "MEWS_score",
              "UVA_score", 
              "shock", 
              "srd", 
              "sevanemia",
              "hivrdtresult", 
              "hivstage34", 
              "newhivdx", 
              "artprior",
              "tmpsmxprior", 
              "malariardtresult", 
              "microtbdx", 
              "urinetblamresult",
              "influenzapcrresult", 
              "hospdeathtransf", 
              "daysspentinhospital",
              "kps70less", 
              "death30d",
              "srsq7",
              "trclust")

#Define categorical variables
catVars <- c("sex", 
             "historyoffever", 
             "nightsweats", 
             "headache", 
             "cough", 
             "diarrhea", 
             "shortnessofbreath",
             "painwithurination", 
             "abxmalarialprior", 
             "temp38plus", 
             "templess36", 
             "ams", 
             "qsofa2p", 
             "qsofa1p", 
             "sirs",
             "shock", 
             "srd", 
             "sevanemia",
             "hivrdtresult", 
             "hivstage34", 
             "newhivdx", 
             "artprior",
             "tmpsmxprior", 
             "malariardtresult", 
             "microtbdx", 
             "urinetblamresult",
             "influenzapcrresult", 
             "hospdeathtransf", 
             "kps70less", 
             "death30d",
             "trclust")

tablesrs7 <- CreateTableOne(vars = listVars, 
                            data = reservesrs4, 
                            factorVars = catVars,
                            strata = "srs7",
                            includeNA = TRUE,
                            addOverall = TRUE)
tablesrs7
summary(tablesrs7)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tablesrs7 <- print(tablesrs7, nonnormal = c("age",
                                            "illnessduration",
                                            "heartrate3", 
                                            "resprate3", 
                                            "sbp3",
                                            "o2sat3",
                                            "qsofa_score",
                                            "sirs",
                                            "MEWS_score",
                                            "UVA_score",
                                            "daysspentinhospital",
                                            "srsq7"),
                   exact = c("sex", 
                             "historyoffever", 
                             "nightsweats", 
                             "headache", 
                             "cough", 
                             "diarrhea", 
                             "shortnessofbreath",
                             "painwithurination", 
                             "abxmalarialprior", 
                             "temp38plus", 
                             "templess36", 
                             "ams", 
                             "qsofa2p", 
                             "qsofa1p", 
                             "sirs",
                             "shock", 
                             "srd", 
                             "sevanemia",
                             "hivrdtresult", 
                             "hivstage34", 
                             "newhivdx", 
                             "artprior",
                             "tmpsmxprior", 
                             "malariardtresult", 
                             "microtbdx", 
                             "urinetblamresult",
                             "influenzapcrresult", 
                             "hospdeathtransf", 
                             "kps70less", 
                             "death30d",
                             "trclust"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tablesrs7), file = 'srs7groups.html')

#Barplot for pathogen and organ dysfunction distribution by 19-gene SRS assignments
library(gmodels)
library(ggplot2)
library(ggthemes)
library(reshape2)

#Dataframe for pathogen distribution 
srsmicro <- data.frame(
  SRS = rep(c("SRS-1", "SRS-2"), 4),
  Micro = c(rep("HIV", 2), rep("TB", 2), rep("Malaria", 2), rep("Influenza", 2)),
  Proportion_of_Patients = c(80,47,40,9,11,26,0,3)
)
srsmicro

#Ok now lets melt the data 
srsmicro.m <- melt(srsmicro)
srsmicro.m

srsmicrobarplot<-ggplot(srsmicro.m, aes(x=factor(Micro, level=c("HIV", "Malaria", "TB", "Influenza")), y=value, fill = SRS)) + 
  geom_bar(stat="identity", colour = "black", position = "dodge")+ 
  scale_fill_manual("SRS", values = c("SRS-1" = "#374e55", "SRS-2" = "#df8f44"))+ 
  labs(y= "Proportion of Patients", x = "")
srsmicrobarplot

srsmicrobarplot <- srsmicrobarplot + theme_bw() + 
  theme(text = element_text(size = 17)) + 
  theme(axis.title = element_text(size = 15), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black")) +  
  theme(legend.title= element_blank())        
srsmicrobarplot

#Dataframe for organ dysfunction distribution 
srsorgan <- data.frame(
  SRS = rep(c("SRS-1", "SRS-2"), 4),
  Organ = c(rep("Shock", 2), rep("Acute resp.\n failure", 2), rep("Encephalopathy", 2), rep("Death at\n 30-days", 2)),
  Proportion_of_Patients = c(23,11,37,20,27,11,48,20)
)
srsorgan

#Ok now lets melt the data 
srsorgan.m <- melt(srsorgan)
srsorgan.m

srsorganbarplot<-ggplot(srsorgan.m, aes(x=factor(Organ, level=c("Shock", "Acute resp.\n failure", "Encephalopathy", "Death at\n 30-days")), y=value, fill = SRS)) + 
  geom_bar(stat="identity", colour = "black", position = "dodge")+ 
  scale_fill_manual("SRS", values = c("SRS-1" = "#374e55", "SRS-2" = "#df8f44"))+ 
  labs(y= "Proportion of Patients", x = "")
srsorganbarplot

srsorganbarplot <- srsorganbarplot + theme_bw() + 
  theme(text = element_text(size = 17, colour = "black")) + 
  theme(axis.title = element_text(size = 15), axis.text.y = element_text(colour = "black"), axis.text.x = element_text(colour = "black")) + 
  theme(legend.title= element_blank())        
srsorganbarplot

#Combine plots
library(ggpubr)
srsplotscombined <- ggarrange(srsmicrobarplot, srsorganbarplot, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
srsplotscombined