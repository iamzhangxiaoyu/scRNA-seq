library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
load('aggr_filter_marker.output.Rdata')

aggr <- NormalizeData(aggr, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(aggr)
aggr <- ScaleData(aggr, features = all.genes)
aggr <- FindVariableFeatures(aggr, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(aggr))
seurat_v3_2000_var_genes=VariableFeatures(aggr)
save(seurat_v3_2000_var_genes,file = 'seurat_v3_2000_var_genes.Rdata')
aggr <- RunPCA(aggr, features = VariableFeatures(object = aggr), verbose = FALSE)
aggr <- JackStraw(aggr, num.replicate = 100, dims = 50)
aggr <- ScoreJackStraw(aggr, dims = 1:50)
JackStrawPlot(aggr, dims = 1:50)
ElbowPlot(aggr, ndims = 50)
aggr <- FindNeighbors(aggr, dims = 1:30)
aggr <- FindClusters(aggr, resolution = 0.2)
aggr <- RunTSNE(aggr, dims = 1:30)
DimPlot(aggr, reduction = "tsne")
save(aggr, file = "aggr for slingshot.Rdata")

source('functions.R')

######aggr
load('aggr_cut_slingshot.Rdata')
aggr <- as.data.frame(aggr@assays$RNA@counts)
head(colnames(aggr))
save(aggr, file = "aggr.Rdata")

####aggr_stages
aggr_stages <- aggr$TimePoints
names(aggr_stages) <- colnames(aggr)
table(aggr_stages)
(dim(aggr))
aggr <- aggr[rowSums(aggr)>0,]
(dim(aggr))
save(aggr_stages, file = "aggr_stages.Rdata")

####aggr_clustering
aggr_clustering <- aggr@meta.data$seurat_clusters
aggr_clustering <- paste("Cluster ", aggr_clustering, sep="")
names(aggr_clustering) <- colnames(aggrs)
table(aggr_clustering)
aggr_clustering[aggr_clustering=="Cluster 1"] <- "Cluster 22"
aggr_clustering[aggr_clustering=="Cluster 2"] <- "Cluster 11"
aggr_clustering[aggr_clustering=="Cluster 22"] <- "Cluster 2"
aggr_clustering[aggr_clustering=="Cluster 11"] <- "Cluster 1"
names(aggr_clustering) <- colnames(aggrs)
table(aggr_clustering)
write.csv(aggr_clustering, file="aggr_clustering.csv")

#####aggrs_data
aggrs_data <- as.matrix(aggrs[rownames(aggrs) %in% seurat_v3_2000_var_genes,])
save(aggrs_data,file = 'aggrs_high_sv_matrix.Rdata')


####DiffusionMap
source('functions.R')
source('colorPalette.R')
load('aggrs_rpkm.Rdata')
load('aggr_stages.Rdata')
load('aggrs_high_sv_matrix.Rdata')
tmp=read.csv('aggr_clustering.csv') 
aggr_clustering=tmp[,2];names(aggr_clustering)=tmp[,1]
table(aggr_clustering)

aggr_dm <- run_diffMap(
  aggrs_data, 
  aggr_clustering,
  sigma=15
)


