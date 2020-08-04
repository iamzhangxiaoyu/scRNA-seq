library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
aggr.data <- Read10X(data.dir = "aggr")
aggr <- CreateSeuratObject(counts = aggr.data, project = "aggr", min.cells = 3, min.features = 200)
timePoints <- sapply(colnames(aggr), function(x) unlist(strsplit(x, "\\-"))[2]) 
timePoints <-ifelse(timePoints == '1', 'Day_2', 
                    ifelse(timePoints == '2', 'Day_4',
                           ifelse(timePoints == '3', 'Day_3', 'Day_0')))
table(timePoints)
aggr <- AddMetaData(object = aggr, metadata = timePoints, col.name = 'TimePoints')
aggr[["percent.mt"]] <- PercentageFeatureSet(aggr, pattern = "^mt-")
VlnPlot(aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'TimePoints')

aggr1 <- subset(aggr, subset = nFeature_RNA > 3000 & nFeature_RNA < 8000 & percent.mt < 7.5 & percent.mt > 5)
VlnPlot(aggr1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'TimePoints')
expr_matrix <- as.matrix(aggr@assays$RNA@counts)

sample_sheet <- data.frame(cells=colnames(aggr@assays$RNA@counts),  
                           cellType=aggr@meta.data)
rownames(sample_sheet)<- colnames(aggr@assays$RNA@counts)
gene_annotation <- as.data.frame(rownames(aggr@assays$RNA@counts))
rownames(gene_annotation)<- rownames(aggr@assays$RNA@counts)
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
HSMM

HSMM <- detectGenes(HSMM, min_expr = 1)
HSMM <- HSMM[fData(HSMM)$num_cells_expressed > 100, ]
HSMM

pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(HSMM), color = cellType.TimePoints, geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < 3e4]
table(HSMM@phenoData@data$cellType.TimePoints)


HSMM <- estimateSizeFactors(HSMM)
#save(HSMM,file = 'HSMM.output.Rdata')
load('HSMM.output.Rdata')
HSMM <- estimateDispersions(HSMM)
HSMM

cds=HSMM

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 5) 

plot_cell_clusters(cds, 1, 2, color = "cellType.TimePoints")

table(pData(cds)$Cluster)

plot_cell_clusters(cds, 1, 2 )

diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Cluster")

sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)

sig_genes$gene_short_name = rownames(sig_genes)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds) 
# 然后降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Cluster")
plot_cell_trajectory(cds, color_by = "cellType.TimePoints")
plot_genes_in_pseudotime(cds[head(sig_genes$gene_short_name),], 
                         color_by = "Cluster")

save(cds, file = 'pseudotime_reduceDimension_no cut cells.Rdata')