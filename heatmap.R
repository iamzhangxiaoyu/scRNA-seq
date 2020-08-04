rm(list=ls())
options(stringsAsFactors = F)

## 加载一系列自定义函数及必备R包，及配色
source('functions_heatmap.R')
# 加载转录组上游数据分析得到的表达矩阵
load('aggr_filter_marker.output.Rdata')
#aggr.markers <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#save(aggr.markers, file = 'aggr.markers.Rdata')
#write.csv(aggr.markers, quote = FALSE, file= "marker.csv")


#top20 <- aggr.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#write.csv(top20, quote = FALSE, file= "top20.csv")

#clustering <- aggr@meta.data$seurat_clusters
#clustering <- paste("Cluster ", clustering, sep="")
#names(clustering) <- colnames(aggr@assays$RNA@counts)
#table(clustering)
#write.csv(clustering, file="clustering.csv")
#clustering[clustering=="Cluster 1"] <- "Cluster 11"
#clustering[clustering=="Cluster 2"] <- "Cluster 22"
#clustering[clustering=="Cluster 11"] <- "Cluster 2"
#clustering[clustering=="Cluster 22"] <- "Cluster 1"
#table(clustering)
#write.csv(clustering, file="clustering.csv")





tmp=read.csv('clustering.csv') 
clustering=tmp[,2];names(clustering)=tmp[,1]
table(clustering)
# 从文章里面拿到的基因列表。
gs=read.csv('marker copy.csv',header = FALSE)[,1]
markerGenes <- gs

# 首先把表达矩阵拆分，重新整理，方便绘图。
gene_subset <- as.matrix(log(aggr@assays$RNA@counts[rownames(aggr@assays$RNA@counts) %in% markerGenes,]+1))
gene_subset[1:4,1:4]
cl0_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="Cluster0"])]
cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="Cluster1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="Cluster2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="Cluster3"])]
#cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="Cluster4"])]

heatmap_gene_subset <- cbind(
  cl0_gene_subset,
  cl1_gene_subset, 
  cl2_gene_subset,
  cl3_gene_subset)
  #cl4_gene_subset


heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- aggr@meta.data$TimePoints
table(heatmap_female_stages)
# 整理好的 heatmap_gene_subset 矩阵 后续绘制热图，共6个发育时期。
# 分成4组。


colbreaks <- c(ncol(cl0_gene_subset),
               ncol(cl0_gene_subset)+ncol(cl1_gene_subset),
               ncol(cl0_gene_subset)+ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
               ncol(cl0_gene_subset)+ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset))
               #ncol(cl0_gene_subset)+ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)+ncol(cl4_gene_subset))

cluster_color <- c(
  Cluster0="#560047",
  Cluster1="#a53bad", 
  Cluster2="#eb6bac", 
  Cluster3="#ffa8a0")
  #Cluster4="#f7dad7"

stage_color=c(
  Day_0="#2754b5", 
  Day_2="#8a00b0", 
  Day_3="#d20e0f", 
  Day_4="#f77f05"
)

tiff(file="DEG_seurat.tiff", 
     res = 300, height = 12, width = 20, units = 'cm')
library(pheatmap)
#tiff('DEG_seurat.tiff', res = 300, width=12, height=8)
#pdf('DEG_seurat.pdf', width=12, height=8)
plot_heatmap_2(
  heatmap_gene_subset, 
  clustering, 
  heatmap_female_stages, 
  rowbreaks, 
  colbreaks,
  cluster_color,
  stage_color
)
dev.off()



