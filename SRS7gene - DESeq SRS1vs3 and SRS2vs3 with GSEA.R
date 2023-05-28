#BiocManager::install(organism, character.only = TRUE)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(DOSE)


# SET THE DESIRED ORGANISM
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#read data 
rawCounts = read.csv("SEP_counts_raw_adult_noPreg_v2.csv", row.names = 1 )
coldata = read.csv("reservesrs4_srs7_assignments.csv")
countdata = rawCounts[, colnames(rawCounts) %in% coldata$pid]


##############################
####FIRST ANALYSIS############
##############################
#SRS 1 vs. SRS3 (SRS3 as reference group)

###transform categories to factors#######################
coldata1 = coldata
countdata1 = countdata
keepC <- coldata1$srs7 %in%  c('SRS1', 'SRS3')
coldata1 <- coldata1[keepC, ]
coldata1$condition <- factor(coldata1$srs7, levels = c("SRS3", "SRS1"))
countdata1 = countdata1[,  colnames(countdata1) %in% coldata1$pid]

dds1 <- DESeqDataSetFromMatrix(countData = countdata1,
                              colData = coldata1,
                              design= ~ condition )
keep <- rowSums(counts(dds1)) >= 10
dds1 <- dds1[keep,]

dds1 <- DESeq(dds1)
resultsNames(dds1) # lists the coefficients
res1 <- results(dds1, name="condition_SRS1_vs_SRS3")
# or to shrink log fold changes association with condition:
res1 <- lfcShrink(dds1, coef="condition_SRS1_vs_SRS3", type="apeglm")
#write to file
out1= data.frame(res1)



#####GSEA
# reading in data from deseq2
df = out1
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- rownames(df)
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse1 <- gseGO(geneList=gene_list, 
             keyType = "SYMBOL", 
             pvalueCutoff = 1,
             eps = 0,
             verbose = TRUE, 
             OrgDb = organism, 
             nPermSimple = 100000
             )

write.csv(as.data.frame(gse1), "SRS7-gsea-SRS1vs3.csv", row.names = FALSE)


##############################
####SECOND ANALYSIS###########
##############################
#SRS2 vs. SRS3 (SRS3 as reference group)

###transform categories to factors#######################
coldata2 = coldata
countdata2 = countdata
keepC2 <- coldata2$srs7 %in%  c('SRS2', 'SRS3')
coldata2 <- coldata2[keepC2, ]
coldata2$condition <- factor(coldata2$srs7, levels = c("SRS3", "SRS2"))
countdata2 = countdata2[,  colnames(countdata2) %in% coldata2$pid]

dds2 <- DESeqDataSetFromMatrix(countData = countdata2,
                               colData = coldata2,
                               design= ~ condition )
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

dds2 <- DESeq(dds2)
resultsNames(dds2) # lists the coefficients
res2 <- results(dds2, name="condition_SRS2_vs_SRS3")
# or to shrink log fold changes association with condition:
res2
res2 <- lfcShrink(dds2, coef="condition_SRS2_vs_SRS3", type="apeglm")
res2

out2= data.frame(res2)

#####GSEA
# reading in data from deseq2
df = out2
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- rownames(df)
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse2 <- gseGO(geneList=gene_list, 
              keyType = "SYMBOL", 
              pvalueCutoff = 1,
              eps = 0,
              verbose = TRUE, 
              OrgDb = organism, 
              nPermSimple = 100000
)

write.csv(as.data.frame(gse2), "SRS7-gsea-SRS2vs3.csv", row.names = FALSE)
