rm(list=ls())
options(stringsAsFactors = F)

## 加载一系列自定义函数及必备R包，及配色
source('functions.R')
# 加载转录组上游数据分析得到的表达矩阵
library(Seurat)
library(monocle)
load(file = 'aggr_filter_marker.output.Rdata')

tmp=read.csv('clustering.csv') 
clustering=tmp[,2];names(clustering)=tmp[,1]  
table(clustering)
count_matrix <- aggr@assays$RNA@counts
stages <- aggr@meta.data$TimePoints
table(stages)

expr_matrix <- as.matrix(count_matrix)
sample_sheet <- data.frame(cells=colnames(count_matrix), 
                           stages=stages, 
                           cellType=clustering)
rownames(sample_sheet)<- colnames(count_matrix)
gene_annotation <- as.data.frame(rownames(count_matrix))
rownames(gene_annotation)<- rownames(count_matrix)
colnames(gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# Create a CellDataSet from the relative expression levels
HSMM <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size()
)

HSMM <- detectGenes(HSMM, min_expr = 5)
# HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 5, ]
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 10, ]

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

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
write.csv(de_genes_table, quote = FALSE, file= "DE_genes_per_clusters_4_groups.csv")
