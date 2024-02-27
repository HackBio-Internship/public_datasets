#============== 001 ============
#install libraries
## TCGA Biolinks
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

##edgeR
BiocManager::install("edgeR")
##limma
BiocManager::install("limma")
##EDASeq
BiocManager::install("EDASeq")
#gplots
install.packages('gplots')



##enable libraries
library("TCGAbiolinks")
library(SummarizedExperiment) #installed alongside TCGABiolinks
library(biomaRt)
library("limma")
library("edgeR")
library('EDASeq')
library(gplots)

#project information
getProjectSummary("TCGA-GBM")
?GDCquery

##Download and preprocess data
gbmQ <- GDCquery(project = "TCGA-GBM",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

GDCdownload(gbmQ)
gbm.data <- GDCprepare(gbmQ)
head(gbm.data)

#explore some metadata information
gbm.data$race
gbm.data$tumor_descriptor
gbm.data$barcode

#we can create a simple metadata for our use case
simpleMeta <- data.frame("barcode" = gbm.data$barcode,
                         "race" = gbm.data$race,
                         'tumor_type' = gbm.data$tumor_descriptor)


























#============== 002 ============
##select unstranded dataset
gbm.raw.data <- assays(gbm.data)
dim(gbm.raw.data$unstranded) #data size

#let's downsize our data to 5 primary and 5 recurring data alone.
selectedBarcodes <- c(subset(simpleMeta, tumor_type == "Recurrence")$barcode[c(1:5)],
                          subset(simpleMeta, tumor_type == "Primary")$barcode[c(1:5)])

selectedData <- gbm.raw.data$unstranded[, c(selectedBarcodes)]
dim(selectedData)

#Data  normalization and filtering
normData <- TCGAanalyze_Normalization(tabDF = selectedData, geneInfo = geneInfoHT, method = "geneLength")
#then filter
filtData <- TCGAanalyze_Filtering(tabDF = normData,
                                  method = "quantile", 
                                  qnt.cut =  0.25)












#============== 003 ============
#we can decide to visualize it as a heatmap, but 60000 is too much and a waste of precious time
#instead, we will pick the top n (50, 100, etc) differentially expressed genes  

#Differential expressed analysis (dea)
selectResults <- 
  TCGAanalyze_DEA(mat1 = filtData[, c(selectedBarcodes)[1:5]],
                mat2 = filtData[, c(selectedBarcodes)[6:10]],
                Cond1type = "Recurrence",
                Cond2type = "Primary",
                pipeline = "edgeR", #can be "limma"
                fdr.cut = 0.01, 
                logFC.cut = 2)

plot(selectResults$logFC, -log10(selectResults$FDR))

#Differential expression analysis with Treatment levels
selectResults.Level <- 
  TCGAanalyze_LevelTab(selectResults,"Recurrence","Primary",
                       filtData[, c(selectedBarcodes)[1:5]],
                       filtData[, c(selectedBarcodes)[6:10]])
  











#============== 004 ============
#now we can visualize with a heatmap
head(selectResults.Level)
dim(selectResults.Level)

heat.data <- filtData[rownames(selectResults.Level),]

#color the plot by the kind of tumor 

cancer.type <- c(rep("Recurrence", 5), rep("Primary", 5))

ccodes <- c()

for (i in cancer.type) {
  if (i == "Recurrence") {
    ccodes <- c(ccodes, "red")
  }else{
    ccodes <- c(ccodes, "blue")
  }
}


heatmap.2(x = as.matrix(heat.data),
          col = hcl.colors(10, palette = 'Blue-Red 2'), #search hcl colors in r
          Rowv = F, Colv = T,
          scale = 'row', 
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Recurrence vs Primary",
          na.color = 'black',
          ColSideColors = ccodes)











#============== 005 ============
#Principal Component Analysis
pca <- TCGAvisualize_PCA(filtData,
                         selectResults.Level, 
                         ntopgenes = 100, 
                         selectedBarcodes[1:5], 
                         selectedBarcodes[6:10])


















#============== 006 ============
#Enrichment Analysis
#View the volcano plot first
plot(x = selectResults.Level$logFC, y = -log10(selectResults.Level$FDR))

upreg.genes <- rownames(subset(selectResults.Level, logFC > 2))
dnreg.genes <- rownames(subset(selectResults.Level, logFC < -2))


#convert  ensemble IDs to gene IDs using biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                      filters = 'ensembl_gene_id', 
                      values = upreg.genes,
                      mart = mart)$hgnc_symbol

dnreg.genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                      filters = 'ensembl_gene_id', 
                      values = dnreg.genes,
                      mart = mart)$hgnc_symbol

#now we can perform enrichment analysis for both
up.EA <- TCGAanalyze_EAcomplete(TFname="Upregulated", upreg.genes)
dn.EA <- TCGAanalyze_EAcomplete(TFname="Downregulated", dnreg.genes)

#now we can visualize our results
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP), 
                        GOBPTab = up.EA$ResBP,
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat,
                        nRGTab = upreg.genes, 
                        nBar = 5,
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)























#============== 007 ============

