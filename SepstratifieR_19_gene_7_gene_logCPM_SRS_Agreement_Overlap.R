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

#Characteristics of patients who switch SRS assignment across SRS19 <--> SRS7
#Create variable for patients who switch assignments
reservesrs4 <- reservesrs4 %>% 
  mutate(switchsrs = if_else(srs19 == srs7, "0", "1"))

#Cross tabs for assignments across 19- and 7-gene classifiers
library(gmodels)
CrossTable(reservesrs4$srs19, reservesrs4$srs7, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)

CrossTable(reservesrs4$switchsrs, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)

#Probabilities of SRS assignment based on 19- or 7-gene classifiers
#19-gene classifier 
reservesrs5 <- merge(reservesrs4, predictions19@SRS_probs, by='row.names',all=TRUE)
rownames(reservesrs5) = reservesrs5$Row.names
reservesrs5$Row.names = NULL
colnames(reservesrs5)[294] = "SRS_1_prob_19"
colnames(reservesrs5)[295] = "SRS_2_prob_19"
colnames(reservesrs5)[296] = "SRS_3_prob_19"

library(dplyr)
reservesrs5 %>% 
  group_by(srs19) %>% 
  summarise(across(SRS_1_prob_19, list(min=min, Q1=~quantile(., probs = 0.25),
                                      median=median, Q3=~quantile(., probs = 0.75),
                                      max=max),  .names = "{.fn}"))

reservesrs5 %>% 
  group_by(srs19) %>% 
  summarise(across(SRS_2_prob_19, list(min=min, Q1=~quantile(., probs = 0.25),
                                       median=median, Q3=~quantile(., probs = 0.75),
                                       max=max),  .names = "{.fn}"))

reservesrs5 %>% 
  group_by(srs19) %>% 
  summarise(across(SRS_3_prob_19, list(min=min, Q1=~quantile(., probs = 0.25),
                                       median=median, Q3=~quantile(., probs = 0.75),
                                       max=max),  .names = "{.fn}"))

#7-gene classifier
reservesrs6 <- merge(reservesrs4, predictions7@SRS_probs, by='row.names',all=TRUE)
rownames(reservesrs6) = reservesrs6$Row.names
reservesrs6$Row.names = NULL
colnames(reservesrs6)[294] = "SRS_1_prob_7"
colnames(reservesrs6)[295] = "SRS_2_prob_7"
colnames(reservesrs6)[296] = "SRS_3_prob_7"

reservesrs6 %>% 
  group_by(srs7) %>% 
  summarise(across(SRS_1_prob_7, list(min=min, Q1=~quantile(., probs = 0.25),
                                       median=median, Q3=~quantile(., probs = 0.75),
                                       max=max),  .names = "{.fn}"))

reservesrs6 %>% 
  group_by(srs7) %>% 
  summarise(across(SRS_2_prob_7, list(min=min, Q1=~quantile(., probs = 0.25),
                                      median=median, Q3=~quantile(., probs = 0.75),
                                      max=max),  .names = "{.fn}"))

reservesrs6 %>% 
  group_by(srs7) %>% 
  summarise(across(SRS_3_prob_7, list(min=min, Q1=~quantile(., probs = 0.25),
                                      median=median, Q3=~quantile(., probs = 0.75),
                                      max=max),  .names = "{.fn}"))

#Cohen kappa for agreement between SRS19 and SRS7 assignments
library(vcd)

assignments <- as.table(rbind(
  c(10,19,1), c(0, 41, 51), c(0, 0, 6)))

categories <- c("SRS1", "SRS2", "SRS3")

dimnames(assignments) <- list(SRS19 = categories, SRS7 = categories)
assignments

res.k <- Kappa(assignments)
res.k

confint(res.k)

#Table of characteristics
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
              "srs19",
              "srs7",
              "srsq19",
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
             "srs19",
             "srs7",
             "trclust")

tableswitchsrs <- CreateTableOne(vars = listVars, 
                                 data = reservesrs4, 
                                 factorVars = catVars,
                                 strata = "switchsrs",
                                 includeNA = TRUE,
                                 addOverall = TRUE)
tableswitchsrs
summary(tableswitchsrs)

#Now specify we need medians(IQR) for continuous variables and Fisher exact for small cell counts
tableswitchsrs <- print(tableswitchsrs, nonnormal = c("age",
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
                                                      "srsq19",
                                                      "srsq7"),
                        exact = c("influenzapcrresult"))

#Now export table to HTML
library(tableHTML)
write_tableHTML(tableHTML(tableswitchsrs), file = 'switchsrs.html')

#Alluvial plot showing the flow of patients across SRS
library(ggplot2)
library(ggalluvial)
library(dplyr)

allusrs <- reservesrs4 %>% 
  group_by(srs19, srs7, switchsrs) %>% 
  summarise(Freq = n())  
allusrs <- na.omit(allusrs)
allusrs$srs19 <- factor(allusrs$srs19)
allusrs$srs7 <- factor(allusrs$srs7)
allusrs$switchsrs <- factor(allusrs$switchsrs)
allusrs

#Set colors
my_colorsallusrs <- RColorBrewer::brewer.pal(9, "Set1")[c(9,4)]

#Make plot
alluvialsrs <- ggplot(allusrs,
                      aes(y = Freq, axis1 = srs19, axis2 = srs7, label = Freq)) +
  geom_alluvium(aes(fill = as.factor(switchsrs)), width = 1/12) +
  geom_stratum(width = 1/4, fill = "grey95", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("SRS-19", "SRS-7"), expand = c(0.05, 0.25)) +
  scale_fill_manual(values = my_colorsallusrs) 

alluvialsrs +  theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(vjust = 6))
alluvialsrs

alluvialsrs <- alluvialsrs + theme_bw() + theme(axis.title.y=element_blank(),
                                                axis.text.y=element_blank(),
                                                axis.ticks.y=element_blank(),
                                                axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(vjust = 6))
alluvialsrs <- alluvialsrs + theme(title = element_text(size = 12),
                                   legend.text = element_text(size = 12),
                                   axis.text.y = element_blank(),
                                   axis.text.x = element_text(size=12))
alluvialsrs

#Now for alluvial plot to show flow of patients who died at 30d
#Lets prepare the data, excluding missing values 
alludeath <- reservesrs4 %>% 
  group_by(death30d, srs19, srs7) %>% 
  summarise(Freq = n()) 
alludeath <- na.omit(alludeath)
alludeath <- alludeath %>% mutate(pct = Freq / sum(Freq))
alludeath$death30d <- factor(alludeath$death30d)
levels(alludeath$death30d) <- c("Alive","Dead")
alludeath

#Set colors
my_colorsallusrs <- RColorBrewer::brewer.pal(9, "Set1")[c(9,3)]

#Plot
alluvialdeath <- ggplot(alludeath,
                      aes(y = Freq, axis1 = death30d, axis2 = srs19, axis3 = srs7, label = Freq)) +
  geom_alluvium(aes(fill = as.factor(death30d)), width = 1/12) +
  geom_stratum(width = 1/5, fill = "grey95", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Vital Status", "SRS-19", "SRS-7"), expand = c(.05, .05)) +
  scale_fill_manual(values = my_colorsallusrs)  

alluvialdeath +  theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(vjust = 3))

alluvialdeath <- alluvialdeath + theme_bw() + theme(axis.title.y=element_blank(),
                                                axis.text.y=element_blank(),
                                                axis.ticks.y=element_blank(),
                                                axis.ticks.x=element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(vjust = 6))
alluvialdeath <- alluvialdeath + theme(title = element_text(size = 12),
                                   legend.text = element_text(size = 12),
                                   axis.text.y = element_blank(),
                                   axis.text.x = element_text(size=12))
alluvialdeath

#Cross tabs between patients assigned to SRS vs. Uganda-derived transcriptional subtypes
library(gmodels)
CrossTable(reservesrs4$trclust, reservesrs4$srs19, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)

CrossTable(reservesrs4$trclust, reservesrs4$srs7, prop.t=TRUE, prop.r=TRUE, prop.c=TRUE)

