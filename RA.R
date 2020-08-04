library(Seurat)
library(dplyr)
library(patchwork)
library(monocle)
library(pheatmap)
library(ggplot2)
source('functions_DEG.R')
RA <- Read10X(data.dir = 'RA')
RA <- CreateSeuratObject(counts = RA, 
                           min.cells = 3, 
                           min.features = 200, 
                           project = "RA")
replicate <- sapply(colnames(RA), function(x) unlist(strsplit(x, "\\-"))[2]) 
table(replicate)
replicate <-ifelse(replicate == "1", "DMSO_Rep1", 
                    ifelse(replicate == "2", "RA_Rep1", 
                           ifelse(replicate == "3", "DMSO_Rep2", "RA_Rep2")))
table(replicate)

RA <- AddMetaData(object = RA, metadata = replicate, col.name = 'replicate')
dim(RA)

fivenum(apply(RA@assays$RNA@counts,1,function(x) sum(x>0) ))
boxplot(apply(RA@assays$RNA@counts,1,function(x) sum(x>0) ))
fivenum(apply(RA@assays$RNA@counts,2,function(x) sum(x>0) ))
hist(apply(RA@assays$RNA@counts,2,function(x) sum(x>0) ))

#RA=RA@assays$RNA@counts[apply(RA@assays$RNA@counts,1, function(x) sum(x>1) > floor(ncol(RA@assays$RNA@counts)/50)),]
dim(RA)

RA[["percent.mt"]] <- PercentageFeatureSet(RA, pattern = "^mt-")
VlnPlot(RA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'replicate')
RA <- subset(RA, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 20)

RA <- NormalizeData(RA, normalization.method = "LogNormalize", scale.factor = 10000)
RA <- FindVariableFeatures(RA, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(RA)
RA <- ScaleData(RA, vars.to.regress = c('nCount_RNA','percent.mt','replicate'), model.use = 'linear', use.umi = FALSE, features = all.genes)
RA <- RunPCA(RA, features = VariableFeatures(object = RA), verbose = FALSE)
VizDimLoadings(RA, dims = 1:2, reduction = "pca")
DimHeatmap(RA, dims = 1:15, cells = 500, balanced = TRUE)
RA <- JackStraw(RA, num.replicate = 100, dims = 50)
RA <- ScoreJackStraw(RA, dims = 1:50)
JackStrawPlot(RA, dims = 1:50)
ElbowPlot(RA, ndims = 50)
RA <- FindNeighbors(RA, dims = 1:30)
RA <- FindClusters(RA, resolution = 0.25)
RA <- RunUMAP(RA, dims = 1:30)
RA <- RunTSNE(RA, dims = 1:30)
save(RA, file = "RA_clusters.Rdata")
load("RA_clusters.Rdata")
DimPlot(RA, reduction = "tsne")
ggsave(filename = 'TSNEPlot_clusters.pdf',width = 10, height =8 )
ggsave(filename = 'TSNEPlot_clusters.tiff',width = 10, height =8 )
DimPlot(RA, reduction = "tsne",split.by = 'replicate',group.by = 'replicate')
ggsave(filename = 'TSNEPlot_clusters_replicate.pdf',width = 40, height =8 )
ggsave(filename = 'TSNEPlot_clusters_replicate.tiff',width = 40, height =8 )
VlnPlot(RA, features = "Ddx4", slot = "counts", log = TRUE, pt.size = 0)
ggsave(filename = 'Violin plot_Ddx4.tiff',width = 10, height =8 )

p1<- FeaturePlot(RA, features = "Stra8", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p1
ggsave(filename = 'TSNEPlot_subcluster_Zbtb16.tiff',width = 10, height =8)


RA <- subset(RA, subset = Ddx4>0)
RA <- FindNeighbors(RA, dims = 1:30)
RA <- FindClusters(RA, resolution = 0.1)
RA <- RunUMAP(RA, dims = 1:30)
RA <- RunTSNE(RA, dims = 1:30)
DimPlot(RA, reduction = "tsne", group.by = 'replicate')
DimPlot(RA, reduction = "tsne")
ggsave(filename = 'TSNEPlot_subclusters.pdf',width = 10, height =8 )
ggsave(filename = 'TSNEPlot_subclusters.tiff',width = 10, height =8 )

p1<- FeaturePlot(RA, features = "Stra8", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p1
ggsave(filename = 'TSNEPlot_subcluster_Zbtb16.tiff',width = 10, height =8)

DimPlot(RA, reduction = "tsne", group.by = 'compare')
ggsave(filename = 'TSNEPlot_subclusters_group.tiff',width = 10, height =8 )
ggsave(filename = 'TSNEPlot_subclusters_group.pdf',width = 10, height =8 )



RA.markers <- FindAllMarkers(RA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
save(RA,file = 'RA_Ddx4.Rdata')

compare <- sapply(colnames(RA), function(x) unlist(strsplit(x, "\\-"))[2]) 
table(compare)
compare <-ifelse(compare == "1", "DMSO", 
                   ifelse(compare == "2", "RA", 
                          ifelse(compare == "3", "DMSO", "RA")))
table(compare)

RA <- AddMetaData(object = RA, metadata = compare, col.name = 'compare')
dim(RA)




load('RA_Ddx4.Rdata')

celltype <-RA$compare
clustering <- RA$compare
names(clustering) <- colnames(RA@assays$RNA@counts)
table(clustering)
write.csv(clustering, file="clustering.csv")
tmp=read.csv('clustering.csv') 
clustering=tmp[,2];names(clustering)=tmp[,1]  
table(clustering)

count_matrix <- RA@assays$RNA@counts
expr_matrix <- as.matrix(count_matrix)
sample_sheet <- data.frame(cells=colnames(count_matrix), 
                           cellType=celltype)
rownames(sample_sheet)<- colnames(count_matrix)
gene_annotation <- as.data.frame(rownames(count_matrix))
rownames(gene_annotation)<- rownames(count_matrix)
colnames(gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

HSMM <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size()
)

HSMM <- detectGenes(HSMM, min_expr = 1)
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
HSMM

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM

diff_test_res <- differentialGeneTest(
  HSMM,
  fullModelFormulaStr="~cellType",
  cores = 4
)

sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)

print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))

diff_test_res <- subset(diff_test_res, qval< 0.01)

save(diff_test_res,file='diff_test_res.Rdata')
save(HSMM, file = 'HSMM.Rdata')
load('HSMM.Rdata')
load('diff_test_res.Rdata')

cluster_nb <- unique(clustering)
mean_per_cluster <- vector()
diff_test_res <- diff_test_res[order(rownames(diff_test_res)),]
count_matrix <- count_matrix[order(rownames(count_matrix)),]
count_de_genes <- count_matrix[rownames(count_matrix) %in% diff_test_res$genes,]
print(dim(count_de_genes))
for (clusters in cluster_nb) {
  # print(head(count_de_genes[,
  # 		colnames(count_de_genes) %in% names(clustering[clustering==clusters])
  # 	]))
  mean <- rowMeans(
    as.matrix(count_de_genes[,
                             colnames(count_de_genes) %in% names(clustering[clustering==clusters])
                             ])
  )
  names(mean) <- clusters
  mean_per_cluster <- cbind(
    mean_per_cluster,
    mean
  )
}
colnames(mean_per_cluster) <- cluster_nb
up_reg_cluster <- colnames(mean_per_cluster)[apply(mean_per_cluster,1,which.max)]
de_genes_table <- data.frame(
  diff_test_res,
  mean_per_cluster,
  cluster=up_reg_cluster
)
write.csv(de_genes_table, quote = FALSE, file= "DE_genes_0.01.csv")





cds=HSMM
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F) 
cds <- reduceDimension(cds, max_components = 2, num_dim =8,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 7) 
plot_cell_clusters(cds, 1, 2, color = "cellType")
plot_cell_clusters(cds, 1, 2 )


ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)

sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
sig_genes$gene_short_name = rownames(sig_genes)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )


plot_cell_trajectory(cds, color_by = "cellType",cell_size = 3,cell_name_size = 3)  
plot_genes_in_pseudotime(cds[head(sig_genes$gene_short_name),], 
                         color_by = "Cluster")

save(cds,file = 'cds.output.Rdata')
load('cds.output.Rdata')
ggsave(filename = 'pseudotime.tiff',width = 10, height =8 )

























load('RA_Ddx4.Rdata')
source('functions_heatmap.R')
#clustering <- Testis@meta.data$celltype
#names(clustering) <- colnames(Testis@assays$RNA@counts)
#table(clustering)
#write.csv(clustering, file="clustering.csv")

tmp=read.csv('clustering.csv') 
clustering=tmp[,2];names(clustering)=tmp[,1]
table(clustering)
# 从文章里面拿到的基因列表。
gs=read.csv('RA response.csv',header = FALSE)[,1]
markerGenes <- gs

# 首先把表达矩阵拆分，重新整理，方便绘图。
gene_subset <- as.matrix(log(RA@assays$RNA@counts[rownames(RA@assays$RNA@counts) %in% markerGenes,]+1))
cl0_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="DMSO"])]
cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(clustering[clustering=="RA"])]

heatmap_gene_subset <- cbind(
  cl0_gene_subset,
  cl1_gene_subset) 

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]

colbreaks <- c(ncol(cl0_gene_subset),
               ncol(cl0_gene_subset)+ncol(cl1_gene_subset))
               
cluster_color <- c(DMSO="pink",
                   RA="grey")

pdf(file="DMSO vs RA_response.pdf", 
     height = 4, width = 8)
library(pheatmap)
#tiff('DEG_seurat.tiff', res = 300, width=12, height=8)
#pdf('DEG_seurat.pdf', width=12, height=8)
plot_heatmap_2(
  heatmap_gene_subset, 
  clustering, 
  colbreaks,
  cluster_color
)
dev.off()




load('RA_Ddx4.Rdata')
VlnPlot(RA, features = "Stra8", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Stra8.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Rarb", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Rarb.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Cyp26a1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Cyp26a1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Cyp26b1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Cyp26b1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Gm4969", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Gm4969.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Prdm9", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Prdm9.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Dmc1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Dmc1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Smc1b", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Smc1b.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Stag3", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Stag3.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Etv5", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Etv5.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Pou5f1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Pou5f1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Id4", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Id4.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Sall4", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Sall4.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Gfra1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Gfra1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Zbtb16", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Zbtb16.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Sycp1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Sycp1.tiff',width = 10, height =8 )

VlnPlot(RA, features = "D6Mm5e", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_D6Mm5e.tiff',width = 10, height =8 )

VlnPlot(RA, features = "Pou5f1", slot = "counts", log = TRUE, pt.size = 0,group.by = 'compare')
ggsave(filename = 'Violin plot_Pou5f1.tiff',width = 10, height =8 )


