rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(Seurat)
library(monocle)
library(ROTS)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
de_genes <-read.csv('DE_genes_per_clusters_4_groups.csv')
gene_names <- subset(de_genes, qval<0.05)
gene_names <- gene_names$genes

# Convert gene ID into entrez genes
entrez_genes <- bitr(gene_names, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]

de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,
                             c("genes", "cluster")]
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
  cluster=de_gene_clusters$cluster
)
table(de_gene_clusters$cluster)
list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                               de_gene_clusters$cluster)
# Run full GO enrichment test
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

# Run GO enrichment test and merge terms 
# that are close to each other to remove result redundancy
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)

# Plot both analysis results
pdf('GO_term_DE_cluster.pdf',width = 11,height = 10)
dotplot(formula_res, showCategory=20)
dev.off()
pdf('GO_term_DE_cluster_simplified.pdf',width = 11,height = 10)
dotplot(lineage1_ego, showCategory=20)
dev.off()


# Save results
write.csv(formula_res@compareClusterResult, 
          file="GO_term_DE_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="GO_term_DE_cluster_simplified.csv")


#########################################################
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichKEGG", 
  organism = "mmu",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

pdf('KEGG_DE_cluster.pdf',width = 11,height = 10)
dotplot(formula_res, showCategory=20)
dev.off()
write.csv(formula_res@compareClusterResult, 
          file="KEGG_DE_cluster.csv")

