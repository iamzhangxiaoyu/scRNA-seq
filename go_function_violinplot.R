library(data.table)
gs=read.csv('mouse oxidative phosphorylation.csv')
gs=as.data.table(gs[!duplicated(gs$Symbol), ]) ######去重复基因#######
write.csv(gs,file = "mouse oxidative phosphorylation.csv")

library(dplyr)
library(Matrix)
library(gplots)
library(ggplot2)
library(Seurat)
load(file = 'aggr_filter_marker.output.Rdata')
aggr$seurat_clusters <- factor(aggr$seurat_clusters, levels = c('0', '2', '1', '3', '4'))
aggr@active.ident<- factor(aggr@active.ident, levels = c('0', '2', '1', '3', '4'))
new.cluster.ids <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4")
names(new.cluster.ids) <- levels(aggr)
aggr <- RenameIdents(aggr, new.cluster.ids)

colP<-c('#560047', '#a53bad', '#eb6bac', '#ffa8a0', '#f7dad7')

x=read.table(paste0("mouse oxidative phosphorylation.txt"),stringsAsFactors=FALSE)[,1]
length(x)
chry.genes <- x[which(x %in% rownames(aggr@assays$RNA@data))]
length(chry.genes)
percent.y <- colSums(expm1(aggr@assays$RNA@data[chry.genes, ]))/length(chry.genes)
aggr <- AddMetaData(aggr, percent.y, "percent.mouseoxidativephosphorylation")
VlnPlot(aggr, features = "percent.mouseoxidativephosphorylation", log = FALSE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin percent.mouse oxidative phosphorylation.pdf',width = 10, height =8 )
