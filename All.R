library(Seurat)
library(tidyverse)
library(ggplot2)
library(monocle)
library(data.table)
library(dplyr)
a=read.table('GSE107746_Folliculogenesis_FPKM.log2.txt.gz',header = T ,sep = '\t',row.names = 'gene')
a=a[,1:80]
Stage <- sapply(colnames(a), function(x) unlist(strsplit(x, "\\_"))[1]) 
Stage <-ifelse(Stage ==  'Primordial', 'Primordial_Oocyte', 
              ifelse(Stage ==  'Primary', 'Primary_Oocyte',      
                     ifelse(Stage ==  'Secondary', 'Secondary_Oocyte', 
                            ifelse(Stage ==  'Antral', 'Antral_Oocyte','Preovulatory_Oocyte'))))  
table(Stage)

Cell <- sapply(colnames(a), function(x) unlist(strsplit(x, "\\O"))[1])
Cell <-ifelse(Cell ==  'Primordial', 'Primordial', 'Oocyte')
table(Cell)

a <- CreateSeuratObject(a, 
                        min.cells = 3, min.features = 200, 
                        project = 'a') 
a
a <- AddMetaData(object = a, metadata = Stage, col.name = 'Stage')
a <- AddMetaData(object = a, metadata = Cell, col.name = 'Cell')
head(a@meta.data)





b=read.table('GSE107746_Folliculogenesis_FPKM.log2.txt.gz',header = T ,sep = '\t',row.names = 'gene')
b=b[,81:151]
Stage <- sapply(colnames(b), function(x) unlist(strsplit(x, "\\_"))[1]) 
Stage <-ifelse(Stage ==  'Primordial', 'Primordial_GC', 
               ifelse(Stage ==  'Primary', 'Primary_GC',      
                      ifelse(Stage ==  'Secondary', 'Secondary_GC', 
                             ifelse(Stage ==  'Antral', 'Antral_GC','Preovulatory_GC'))))  
table(Stage)

Cell <- sapply(colnames(b), function(x) unlist(strsplit(x, "\\O"))[1])
Cell <-ifelse(Cell ==  'Primordial', 'Primordial', 'Granulosa')
table(Cell)

b <- CreateSeuratObject(b, 
                        min.cells = 3, min.features = 200, 
                        project = 'b') 
b
b <- AddMetaData(object = b, metadata = Stage, col.name = 'Stage')
b <- AddMetaData(object = b, metadata = Cell, col.name = 'Cell')
head(b@meta.data)

All <- merge(a, y = b, add.cell.ids = c("a", "b"), project = "All")

All <- NormalizeData(All, normalization.method = "LogNormalize", scale.factor = 10000)
All <- NormalizeData(All)
All <- FindVariableFeatures(All, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(All)
All <- ScaleData(All, features = all.genes)
All <- RunPCA(All, features = VariableFeatures(object = All))
All <- JackStraw(All, num.replicate = 100)
All <- ScoreJackStraw(All, dims = 1:50)
All <- FindNeighbors(All, dims = 1:30)
All <- FindClusters(All, resolution = 2)
All <- RunUMAP(All, dims = 1:30)
All <- RunTSNE(All, dims = 1:30)
DimPlot(All, reduction = "umap", label = FALSE,pt.size = 3)
DimPlot(All, reduction = "umap", label = FALSE,pt.size = 3,group.by = 'Stage')
DimPlot(All, reduction = "umap", label = FALSE,pt.size = 3,group.by = 'Cell')
All
p1<- FeaturePlot(All, features = "ULK1", reduction = "umap", pt.size = 3, sort.cell = TRUE, label = FALSE, cols = c("lightgrey", "red"))
p1
save(All, file = "All.Rdata")

#####################Oocyte###################################
Oocyte=read.table('GSE107746_Folliculogenesis_FPKM.log2.txt.gz',header = T ,sep = '\t',row.names = 'gene')
Oocyte=Oocyte[,1:80]
Stage <- sapply(colnames(Oocyte), function(x) unlist(strsplit(x, "\\_"))[1]) 
Stage <-ifelse(Stage ==  'Primordial', 'Primordial_Oocyte', 
               ifelse(Stage ==  'Primary', 'Primary_Oocyte',      
                      ifelse(Stage ==  'Secondary', 'Secondary_Oocyte', 
                             ifelse(Stage ==  'Antral', 'Antral_Oocyte','Preovulatory_Oocyte'))))  
table(Stage)

Cell <- sapply(colnames(Oocyte), function(x) unlist(strsplit(x, "\\O"))[1])
Cell <-ifelse(Cell ==  'Primordial', 'Primordial', 'Oocyte')
table(Cell)

Oocyte <- CreateSeuratObject(Oocyte, 
                        min.cells = 3, min.features = 200, 
                        project = 'Oocyte') 
Oocyte
Oocyte <- AddMetaData(object = Oocyte, metadata = Stage, col.name = 'Stage')
Oocyte <- AddMetaData(object = Oocyte, metadata = Cell, col.name = 'Cell')
head(Oocyte@meta.data)

Oocyte <- NormalizeData(Oocyte, normalization.method = "LogNormalize", scale.factor = 10000)
Oocyte <- NormalizeData(Oocyte)
Oocyte <- FindVariableFeatures(Oocyte, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Oocyte)
Oocyte <- ScaleData(Oocyte, features = all.genes)
Oocyte <- RunPCA(Oocyte, features = VariableFeatures(object = Oocyte))
Oocyte <- JackStraw(Oocyte, num.replicate = 100)
Oocyte <- ScoreJackStraw(Oocyte, dims = 1:20)
JackStrawPlot(Oocyte, dims = 1:15)
ElbowPlot(Oocyte)
Oocyte <- FindNeighbors(Oocyte, dims = 1:20)
Oocyte <- FindClusters(Oocyte, resolution = 1.5)
Oocyte <- RunUMAP(Oocyte, dims = 1:20)
DimPlot(Oocyte, reduction = "umap", label = FALSE,pt.size = 3)
DimPlot(Oocyte, reduction = "umap", label = FALSE,pt.size = 3,group.by = 'Stage')
DimPlot(Oocyte, reduction = "umap", label = FALSE,pt.size = 3,group.by = 'Cell')
Oocyte
save(Oocyte, file = "Oocyte.Rdata")

#####################SCENIC###################################
load("Oocyte.Rdata")
exprMat <- as.matrix(Oocyte@assays$RNA@counts)
dim(exprMat)
cellInfo <- data.frame(cellCluster=Oocyte@meta.data[,1:3])
colnames(cellInfo)[1:3] <- c("cellCluster","nUMI","nGene")

table(cellInfo)
saveRDS(cellInfo, file="cellInfo.Rds")
colVars <- list(cellCluster=c("Primordial"="forestgreen", 
                                "Primary"="darkorange", 
                                "Secondary"="magenta4", 
                                "Antral"="hotpink", 
                                "Preovulatory"="red3"))
colVars$cellCluster <- colVars$cellCluster[intersect(names(colVars$cellCluster), cellInfo$cellCluster)]
saveRDS(colVars, file="colVars.Rds")
plot.new(); legend(0,1, fill=colVars$cellCluster, legend=names(colVars$cellCluster))
library(SCENIC)
org <- "hgnc" # or mgi, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "Oocyte" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"

# Save to use at a later time...
saveRDS(scenicOptions, file="scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)
rm(exprMat)
runCorrelation(exprMat_filtered, scenicOptions)

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

exprMat <- as.matrix(Oocyte@assays$RNA@counts)
logMat <- log2(exprMat+1)
dim(exprMat)
library(SCENIC)
scenicOptions <- readRDS("scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
scenicOptions@settings$defaultTsne$dims <- 26
scenicOptions@settings$defaultTsne$perpl <- 26
# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["500bp"]

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod= NULL)
runSCENIC_3_scoreCells(scenicOptions, logMat)


aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="scenicOptions.Rds")
# scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)


nPcs <- c(5) # For toy dataset
# nPcs <- c(5,15,50)

scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="cellCluster", cex=.5)


# Using only "high-confidence" regulons (normally similar)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="cellCluster", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("OTX2", "JUNB", "CEBPZ","ATF6")],], plots="Expression")

# Save AUC as PDF:
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

#par(bg = "black")
par(mfrow=c(1,2))

regulonNames <- c( "OTX2", "CEBPZ")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(0, 10, attr(cellCol,"red"), col="red", cex=.7, pos=4)

text(-20,-10, attr(cellCol,"green"), col="green3", cex=.7, pos=4)

regulonNames <- list(red=c("OTX2", "CEBPZ"),
                     green=c("JUNB"),
                     blue=c( "ATF6"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")

text(5, 15, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(5, 15-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(5, 15-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)

regulons <- loadInt(scenicOptions, "regulons")
regulons[c("OTX2", "CEBPZ")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="OTX2" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="OTX2"]
viewMotifs(tableSubset) 

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_bycellCluster <- sapply(split(rownames(cellInfo), cellInfo$cellCluster),
                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_bycellCluster_Scaled <- t(scale(t(regulonActivity_bycellCluster), center = T, scale=T))

pheatmap::pheatmap(regulonActivity_bycellCluster_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,filename="regulonActivity_byseuratCluster.pdf", width=10, height=20)

topRegulators <- reshape2::melt(regulonActivity_bycellCluster_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)



minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_bycellCluster_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$cellCluster), 
                                                    function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_bycellCluster_Binarized[which(rowSums(regulonActivity_bycellCluster_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5, 
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,filename="regulonActivityBinary_byCellType.pdf", width=10, height=20)

topRegulators <- reshape2::melt(regulonActivity_bycellCluster_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)

library(Seurat)
dr_coords <- Embeddings(Oocyte, reduction="tsne")

tfs <- c("OTX2", "CEBPZ", "JUNB","ATF6")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
