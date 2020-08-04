load('aggr_slingshot_8k.Rdata')
aggrs <- as.data.frame(aggr_cut@assays$RNA@counts)
#####################################################################
rm(list=ls())
options(stringsAsFactors = F)
source('functions.R')
load('aggrs_rpkm.Rdata')
head(colnames(aggrs))
aggr_stages <- sapply(strsplit(colnames(aggrs), "_"), `[`, 1)
library(stringr)
aggr_stages=str_split(colnames(aggrs), "_",simplify = T)[,1]
names(aggr_stages) <- colnames(aggrs)
table(aggr_stages)
(dim(aggrs))
aggrs <- aggrs[rowSums(aggrs)>0,]
(dim(aggrs))
exprSet=aggrs
mean_per_gene <- apply(exprSet , 1, mean, na.rm = TRUE) #对表达矩阵每行求均值
sd_per_gene <- apply(exprSet, 1, sd, na.rm = TRUE) #对表达矩阵每行求标准差
mad_per_gene <-   apply(exprSet, 1, mad, na.rm = TRUE) #对表达矩阵每行求绝对中位差
cv = sd_per_gene/mean_per_gene

library(matrixStats)
var_per_gene <- rowVars(as.matrix(exprSet))
cv2=var_per_gene/mean_per_gene^2
cv_per_gene <- data.frame(mean = mean_per_gene,
                          sd = sd_per_gene,
                          mad=mad_per_gene,
                          var=var_per_gene,
                          cv=cv,
                          cv2=cv2)
rownames(cv_per_gene) <- rownames(exprSet)
head(cv_per_gene)
cv_per_gene=cv_per_gene[cv_per_gene$mean>1,]
with(cv_per_gene,plot(log10(mean),log10(cv2)))

library(psych)
pairs.panels(cv_per_gene, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)


aggrs_data <- getMostVarGenes(aggrs, fitThr=1)
aggrs_data <- log(aggrs_data+1)
dim(aggrs_data)
dim(aggrs_data)
aggrs_data[1:4,1:4]
save(aggrs_data,file = 'aggrs_high_sv_matrix.Rdata')


aggr_sub_pca <- FactoMineR::PCA(
  t(aggrs_data), 
  ncp = ncol(aggrs_data), 
  graph=FALSE
)

# Estimate which PCs contain significant information with Jackstraw
# the result is 9 significant PCs
# 挑选显著的主成分数量，在seurat包有包装好的函数。
significant_pcs <- jackstraw::permutationPA(
  aggr_sub_pca$ind$coord, 
  B = 100, 
  threshold = 0.05, 
  verbose = TRUE, 
  seed = NULL
)$r
significant_pcs


# Compute and plot the t-SNE using the significant PCs
# 基于PCA的结果挑选显著的主成分进行tSNE
# 主要是因为 tSNE 耗费计算量。
aggr_t_sne <- run_plot_tSNE(
  pca=aggr_sub_pca,
  pc=significant_pcs,
  iter=5000,
  conditions=aggr_stages,
  colours=aggr_stagePalette
)
ggsave('tSNE_by_stage.pdf')
save(aggr_t_sne,file='step1-aggr_t_sne.Rdata')

# Recompute PCA with FactomineR package for the clustering
res.pca <- PCA(
  t(aggrs_data), 
  ncp = significant_pcs, 
  graph=FALSE
)

# Clustering cells based on the PCA with 
# a minimum of 4 expected clusters
# Hierarchical Clustering On Principle Components (HCPC)
# 基于PCA的结果进行层次聚类，这里可以指定是4类。
res.hcpc <- HCPC(
  res.pca, 
  graph = FALSE,
  min=4
)

# Plot the hierarchical clustering result
plot(res.hcpc, choice ="tree", cex = 0.6)

# Extract clustering results
aggr_clustering <- res.hcpc$data.clust$clust
aggr_clustering <- paste("C", aggr_clustering, sep="")
names(aggr_clustering) <- rownames(res.hcpc$data.clust)
table(aggr_clustering)
write.csv(aggr_clustering, file="aggr_clustering.csv")
aggr_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)



# Plot the t-SNE colored by cell clusters
head(aggr_t_sne)
# 包装了一个ggplot绘图函数。
aggr_t_sne_new_clusters <- plot_tSNE(
  tsne=aggr_t_sne, 
  conditions=aggr_clustering, 
  colours= aggr_clusterPalette
)
ggsave('tSNE_cluster.pdf')


rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(densityClust)
library(destiny)
library(rgl)
source('functions.R')
source('colorPalette.R')

# 载入RPKM形式的表达矩阵。
load('aggrs_rpkm.Rdata')
# Extract embryonic stage data from cell names
head(colnames(aggrs))
aggr_stages <- sapply(strsplit(colnames(aggrs), "_"), `[`, 1)
names(aggr_stages) <- colnames(aggrs)
table(aggr_stages)

### 载入第一步挑选到的800多个高变异基因
load(file = 'aggrs_high_sv_matrix.Rdata')
dim(aggrs_data)

## 载入第一步的tSNE后细胞类群信息
tmp=read.csv('aggr_clustering.csv') 
aggr_clustering=tmp[,2];names(aggr_clustering)=tmp[,1]
names(aggr_clustering)=colnames(aggrs)
table(aggr_clustering)

# 下面运行 DiffusionMap 需要822个基因的表达矩阵
# 以及它们被tSNE+DBSCAN聚的4类

# 可以把 DiffusionMap 类比 PCA 分析

# Compute the Diffusion map
aggr_dm <- run_diffMap(
  aggrs_data, 
  aggr_clustering,
  sigma=15
)

save(aggr_dm,aggr_clustering,aggr_stages,
     file = 'diffusionMap_output.Rdata')
# Plot the Eigen values per diffusion component (similar to screeplots for PCAs)
plot_eigenVal(
  dm=aggr_dm
)

save(aggr_dm,aggr_lineage,aggr_pseudotime,aggrs,aggrs_data,aggr_clustering,aggr_stages, file = 'Pseudotime_for_slingshot.Rdata')
load('Pseudotime_for_slingshot.Rdata')
library(ggplot2)
library(data.table)
library(densityClust)
library(destiny)
library(rgl)
source('functions.R')
source('colorPalette.R')
pseudotime_lin <- aggr_pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime

aggr_pseudotime[,"curve1"] <- pseudotime_lin1_percent

# Twisted palette to plot gene expression by lineages along pseudotime
aggr_clusterPalette2 <- c(
  "#ff6663", 
  "#3b3561"
)

# Plot one gene expression over pseudotime for the lineage 1
plot_smoothed_gene_per_lineage(
  rpkm_matrix=aggrs, 
  pseudotime=aggr_pseudotime, 
  lin=c(1),
  gene="Sycp1", 
  stages=aggr_stages, 
  clusters=aggr_clustering, 
  stage_colors=aggr_stagePalette,
  cluster_colors=aggr_clusterPalette,
  lineage_colors=aggr_clusterPalette2
)
ggsave("Sycp1.pdf",width = 10, height = 6)

gs=read.csv('all.csv',header = FALSE)[,1]
gene_list<- gs[which(gs %in% rownames(aggrs))]

gene_list <- c("Bbx",
               "Arid4a",
              "Arid4b",
              "Atoh8",
               "Bcl11a",
               "Creb3l4",
               "Foxs1",
               "Lhx9",
               "Mxd4",
               "Nr1h3",
               "Pax9",
               "Zfp780b",
               "Zfp867",
               "Zfp157",
               "4930522L14Rik",
               "Pou2f2",
               "Stat1",
               "Tox",
               "Zfp932",
               "Ascl2",
               "Epas1",
               "Klf4",
              "Mga",
               "Nfib",
               "Nfix",
               "Sohlh1",
               "Sox3",
               "Tshz2",
               "Zfp292" )    


plot_smoothed_genes <- function(genes, lin){
  aggr_clusterPalette2 <- c("#ff6663", "#3b3561")
  for (gene in genes){
    plot_smoothed_gene_per_lineage(
      rpkm_matrix=aggrs, 
      pseudotime=aggr_pseudotime, 
      lin=lin,
      gene=gene, 
      stages=aggr_stages, 
      clusters=aggr_clustering, 
      stage_colors=aggr_stagePalette,
      cluster_colors=aggr_clusterPalette,
      lineage_colors=aggr_clusterPalette2
    )
  }
}

pdf("Pseu_tf.pdf", width=10, height=6)
plot_smoothed_genes(gene_list, 1) # plot only lineage 1
#plot_smoothed_genes(gene_list, 2) # plot only lineage 2
#plot_smoothed_genes(gene_list, c(1,2)) # plot the two moleages in the same graph to see the divergence
dev.off()