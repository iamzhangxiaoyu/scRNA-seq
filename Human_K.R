library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)

K1.data <- Read10X(data.dir = "/Users/zhangxiaoyu/Desktop/human kidney/GSE131685_RAW/kidney1/")
K1 <- CreateSeuratObject(counts = K1.data, project = "kidney1", min.cells = 8, min.features = 200)
K2.data <- Read10X(data.dir = "/Users/zhangxiaoyu/Desktop/human kidney/GSE131685_RAW/kidney2/")
K2 <- CreateSeuratObject(counts = K2.data, project = "kidney2", min.cells = 6, min.features = 200)
K3.data <- Read10X(data.dir = "/Users/zhangxiaoyu/Desktop/human kidney/GSE131685_RAW/kidney3/")
K3 <- CreateSeuratObject(counts = K3.data, project = "kidney3", min.cells = 10, min.features = 200)
kid <- merge(x = K1, y = list(K2, K3)) 

kid[["percent.mt"]] <- PercentageFeatureSet(kid, pattern = "^MT-") 
VlnPlot(kid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

kid <- subset(kid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30) 
kid <- NormalizeData(kid, normalization.method = "LogNormalize", scale.factor = 10000)
kid <- NormalizeData(kid) #标准化
kid <- FindVariableFeatures(kid, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(kid), 10)
plot1 <- VariableFeaturePlot(kid)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
kid <- CellCycleScoring(kid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.genes <- rownames(kid)
kid <- ScaleData(kid, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

#Eliminate batch effects with harmony and cell classification
kid <- RunPCA(kid, pc.genes = kid@var.genes, npcs = 20, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
kid <- kid %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)　
harmony_embeddings <- Embeddings(kid, 'harmony')
harmony_embeddings[1:5, 1:5]
kid <- kid %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.25) %>%
  identity()
new.cluster.ids <- c(0,1, 2, 3, 4, 5, 6, 7,8,9,10)
names(new.cluster.ids) <- levels(kid)
kid <- RenameIdents(kid, new.cluster.ids)

#Calculating differentially expressed genes (DEGs) and Save rds file
kid.markers <- FindAllMarkers(kid, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(kid.markers,sep="\t",file="kidney_marker.xls")
save(kid,file="kidney.Rdata")

#Some visual figure generation
DimPlot(kid, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')

DimPlot(kid, reduction = "umap", group.by = "Phase", pt.size = .1)

DimPlot(kid, reduction = "umap", label = TRUE, pt.size = .1)
p1<- FeaturePlot(kid, features = "STRA8", reduction = "umap", pt.size = 3, sort.cell = TRUE, label = FALSE, cols = c("lightgrey", "red"))
p1

p2<- FeaturePlot(kid, features = "SLC22A7", reduction = "umap", pt.size = 3, sort.cell = TRUE, label = FALSE, cols = c("lightgrey", "red"))
p2
