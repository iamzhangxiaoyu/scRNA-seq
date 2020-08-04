rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(Seurat)
library(monocle)
library(ROTS)
library(clusterProfiler)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

###########################################
#                                         #
#             Load files                  #
#                                         #
###########################################
de_genes <-read.csv('marker cluster0.csv')
gene_names <- de_genes$gene

###########################################
#                                         #
#          Transfer ID symbol             #
#                                         #
###########################################
entrez_genes <- bitr(gene_names, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]
#write.csv(entrez_genes, file="entrez_genes.csv")

###########################################
#                                         #
#                 BP                      #
#                                         #
###########################################
ego <- enrichGO(gene = entrez_genes$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.9,
                qvalueCutoff  = 0.9)
ego=DOSE::setReadable(ego, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
head(ego)
pdf('GO_term_BP_marker0.pdf',width = 11,height = 10)
dotplot(ego, showCategory=50)
dev.off()
tiff('GO_term_BP_marker0.tiff',width = 3500,height = 3000, units = "px",res = 300)
dotplot(ego, showCategory=50)
dev.off()
write.csv(ego, file="GO_term_BP_marker0.csv")

###########################################
#                                         #
#                  CC                     # 
#                                         #
###########################################
ego <- enrichGO(gene = entrez_genes$ENTREZID,
                OrgDb = org.Mm.eg.db,
                  ont = "CC", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.9,
                qvalueCutoff  = 0.9)
ego=DOSE::setReadable(ego, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
head(ego)
pdf('GO_term_CC_marker0.pdf',width = 11,height = 10)
dotplot(ego, showCategory=50)
dev.off()
tiff('GO_term_CC_marker0.tiff',width = 3500,height = 3000, units = "px",res = 300)
dotplot(ego, showCategory=50)
dev.off()
write.csv(ego, file="GO_term_CC_marker0.csv")

###########################################
#                                         #
#                  MF                     #
#                                         #
###########################################
ego <- enrichGO(gene = entrez_genes$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "MF", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.9,
                qvalueCutoff  = 0.9)
ego=DOSE::setReadable(ego, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
head(ego)
pdf('GO_term_MF_marker0.pdf',width = 11,height = 10)
dotplot(ego, showCategory=50)
dev.off()
tiff('GO_term_CC_marker0.tiff',width = 3500,height = 3000, units = "px",res = 300)
dotplot(ego, showCategory=50)
dev.off()
write.csv(ego, file="GO_term_MF_marker0.csv")


###########################################
#                                         #
#                  KEGG                   #
#                                         #
###########################################
kk <- enrichKEGG(gene = entrez_genes$ENTREZID,
                   organism = "mmu",
                   keyType = "kegg", 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.9,
                   qvalueCutoff =0.9)
kk=DOSE::setReadable(kk, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
head(kk)
pdf('KEGG_marker0.pdf',width = 11,height = 10)
dotplot(kk, showCategory=50)
dev.off()
tiff('KEGG_marker0.tiff',width = 3500,height = 3000, units = "px",res = 300)
dotplot(kk, showCategory=50)
dev.off()
write.csv(kk, file="KEGG_marker0.csv")
