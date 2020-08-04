library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load(file = 'aggr_filter_marker.output.Rdata')
aggr$seurat_clusters <- factor(aggr$seurat_clusters, levels = c('0', '2', '1', '3', '4'))
aggr@active.ident<- factor(aggr@active.ident, levels = c('0', '2', '1', '3', '4'))
new.cluster.ids <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4")
names(new.cluster.ids) <- levels(aggr)
aggr <- RenameIdents(aggr, new.cluster.ids)
aggr <- subset(aggr, subset = Ddx4 > 0)


subcluster = c("Cluster0","Cluster3")
sub <- subset(aggr, idents = subcluster)
df.data <- GetAssayData(object = sub, slot = "data")

df.group <- data.frame(umi = names(Idents(sub)), 
                       cluster = as.character(sub@meta.data$RNA_snn_res.0.2), 
                       stringsAsFactors = F)
#save(df.data,file = "df.data.Rdata")
#save(df.group,file = "df.group.Rdata")
#load("df.data.Rdata")
#load("df.group.Rdata")
genesets <- getGmt("MousePath_GO_gmt.gmt")
gsvascore <- gsva(data.matrix(df.data), genesets, parallel.sz = 2)
gsvascore[1:5, 1:5]

ha.t <- HeatmapAnnotation(Cluster = df.group$cluster)
Heatmap(as.matrix(gsvascore), 
        show_column_names = F, 
        cluster_rows = T, 
        cluster_columns = T, 
        top_annotation = ha.t, 
        column_split = df.group$cluster, 
        row_names_gp = gpar(fontsize = 8), 
        row_names_max_width = max_text_width(rownames(gsvascore), 
                                             gp = gpar(fontsize = 8)))


### 3. 挑选出要使用小提琴图展示的通路
# 使用limma包进行差异分析，然后根据差异分析结果，对通路进行选择。
cluster.group <- as.numeric(df.group$cluster) 
# 确定分组信息
group2 <- 3 
# 根据子集中细胞类型的名称来定
# 这里cluster6作为实验组（group2)，cluster8作为对照组（group1)
cluster.group[which(cluster.group != group2)] <- 0
cluster.group[which(cluster.group == group2)] <- 1
design <- cbind(sampleGroup1=1, sampleGroup2vs1=cluster.group)

# 利用limma包进行差异分析
fit <- lmFit(gsvascore, design)    
fit <- eBayes(fit)
sigPathways <- topTable(fit, coef="sampleGroup2vs1", 
                        number=Inf, p.value=0.05, adjust="BH")
# 添加一下分组信息，方便导出的数据区分实验组和对照组。
sigPathways$celltype <- rep(group2, length(nrow(sigPathways))) 
# 保存到文件
write.csv(sigPathways, file = "Aggr_DE_genesets_GO_cluster3_vs_0.csv")



# 这里选第一个pathway画图
count <- gsvascore[rownames(sigPathways)[57], , drop = FALSE]
count <- as.data.frame(t(count))
colnames(count) <- "geneset"
count$cluster <- as.character(Idents(sub))
# 用通路的名称作为图形标题
title.name = rownames(sigPathways)[57]

# 添加P值
# 得到cluster6和cluster10中的gsva score 最大值
count.geneset.group1 <- count$geneset[count$cluster == subcluster[1]]
count.geneset.group2 <- count$geneset[count$cluster == subcluster[2]]

# 确定P值添加的位置
ysegment1 <- max(count.geneset.group1)
ysegment2 <- max(count.geneset.group2)
ysegment.max <- max(ysegment1, ysegment2)

# 根据P值确定加几颗星
pval <- sigPathways$P.Value[1]

if (pval < 0.001) {
  pval.label = "***"
} else if (pval < 0.005) {
  pval.label = "**"
} else if (pval < 0.05) {
  pval.laben = "*"
} else if (pval >= 0.05) {
  pval.label = "non.sig"
}

# 自定义颜色
blue <- "#619CD6"
green <- "#89C32E"
p <- ggplot(count, aes(x = cluster, y = geneset, fill = cluster)) +
  geom_violin() +
  scale_fill_manual(values = c(blue, green)) + # 用自定义颜色填充
  theme_classic() +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
        axis.title.x = element_text(color = "black", size = 20), 
        axis.title.y = element_blank(), 
        axis.text = element_text(color = "black", size = 16), 
        axis.line = element_line(colour = "black", size = 0.6), 
        plot.title = element_text(size = 14, hjust = 0.5)) + 
  # 添加图形标题
  ggtitle(title.name) +
  guides(fill = F)
p
p + 
  # 图形中的一个横线和两个竖线
  annotate("segment", x = 1, xend = 2, y = ysegment.max + 0.02, yend = ysegment.max + 0.02) + 
  annotate("segment", x = 1, xend = 1, y = ysegment1 + 0.01, yend = ysegment.max + 0.02) +
  annotate("segment", x = 2, xend = 2, y = ysegment2 + 0.01, yend = ysegment.max + 0.02) +
  # 添加P值对应的星号
  annotate("text", 
           size = 12, # *的大小
           x = 1.5, 
           y = ysegment.max - 0.03, #可以微调*所在的位置
           label = pval.label)

ggsave("PBMC_DE_genesets_hallmarker_cluster6_vs_8.pdf", width = 5.5, height = 5)
