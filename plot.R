library(Seurat)
library(ggplot2)
load(file = 'aggr_filter_marker.output.Rdata')
p1<- FeaturePlot(aggr, features = "Gfra1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p1
ggsave(filename = 'TSNEPlot_aggr_filter_Gfra1.tiff',width = 10, height =8)

p2<- FeaturePlot(aggr, features = "Zbtb16", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p2
ggsave(filename = 'TSNEPlot_aggr_filter_Zbtb16.tiff',width = 10, height =8)

p3<- FeaturePlot(aggr, features = "Nanos3", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p3
ggsave(filename = 'TSNEPlot_aggr_filter_Nanos3.tiff',width = 10, height =8)

p4<- FeaturePlot(aggr, features = "Etv5", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p4
ggsave(filename = 'TSNEPlot_aggr_filter_Etv5.tiff',width = 10, height =8)

p5<- FeaturePlot(aggr, features = "Id4", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p5
ggsave(filename = 'TSNEPlot_aggr_filter_Id4.tiff',width = 10, height =8)

p6<- FeaturePlot(aggr, features = "Kit", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p6
ggsave(filename = 'TSNEPlot_aggr_filter_Kit.tiff',width = 10, height =8)

p7<- FeaturePlot(aggr, features = "Stra8", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p7
ggsave(filename = 'TSNEPlot_aggr_filter_Stra8.tiff',width = 10, height =8)

p8<- FeaturePlot(aggr, features = "Cyp26a1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p8
ggsave(filename = 'TSNEPlot_aggr_filter_Cyp26a1.tiff',width = 10, height =8)

p9<- FeaturePlot(aggr, features = "Rhox13", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p9
ggsave(filename = 'TSNEPlot_aggr_filter_Rhox13.tiff',width = 10, height =8)

p10<- FeaturePlot(aggr, features = "Utf1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p10
ggsave(filename = 'TSNEPlot_aggr_filter_Utf1.tiff',width = 10, height =8)

p11<- FeaturePlot(aggr, features = "Sohlh1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p11
ggsave(filename = 'TSNEPlot_aggr_filter_Sohlh1.tiff',width = 10, height =8)

p12<- FeaturePlot(aggr, features = "Sohlh2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p12
ggsave(filename = 'TSNEPlot_aggr_filter_Sohlh2.tiff',width = 10, height =8)

p13<- FeaturePlot(aggr, features = "Hormad1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p13
ggsave(filename = 'TSNEPlot_aggr_filter_Hormad1.tiff',width = 10, height =8)

p14<- FeaturePlot(aggr, features = "Dmc1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.75)
p14
ggsave(filename = 'TSNEPlot_aggr_filter_Dmc1.tiff',width = 10, height =8)

p15<- FeaturePlot(aggr, features = "Meiob", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p15
ggsave(filename = 'TSNEPlot_aggr_filter_Meiob.tiff',width = 10, height =8)

p16<- FeaturePlot(aggr, features = "Spo11", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p16
ggsave(filename = 'TSNEPlot_aggr_filter_Spo11.tiff',width = 10, height =8)

p17<- FeaturePlot(aggr, features = "Sycp1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p17
ggsave(filename = 'TSNEPlot_aggr_filter_Sycp1.tiff',width = 10, height =8)

p18<- FeaturePlot(aggr, features = "Sycp3", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p18
ggsave(filename = 'TSNEPlot_aggr_filter_Sycp3.tiff',width = 10, height =8)

p19<- FeaturePlot(aggr, features = "Mybl1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p19
ggsave(filename = 'TSNEPlot_aggr_filter_Mybl1.tiff',width = 10, height =8)

p20<- FeaturePlot(aggr, features = "Tex15", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p20
ggsave(filename = 'TSNEPlot_aggr_filter_Tex15.tiff',width = 10, height =8)

p21<- FeaturePlot(aggr, features = "Smc1b", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2.5)
p21
ggsave(filename = 'TSNEPlot_aggr_filter_Smc1b.tiff',width = 10, height =8)

p22<- FeaturePlot(aggr, features = "Meioc", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p22
ggsave(filename = 'TSNEPlot_aggr_filter_Meioc.tiff',width = 10, height =8)

p23<- FeaturePlot(aggr, features = "Tex19.2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p23
ggsave(filename = 'TSNEPlot_aggr_filter_Tex19.2.tiff',width = 10, height =8)

p24<- FeaturePlot(aggr, features = "Tex19.1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p24
ggsave(filename = 'TSNEPlot_aggr_filter_Tex19.1.tiff',width = 10, height =8)

p25<- FeaturePlot(aggr, features = "Tex101", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p25
ggsave(filename = 'TSNEPlot_aggr_filter_Tex101.tiff',width = 10, height =8)

p26<- FeaturePlot(aggr, features = "Spata22", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p26
ggsave(filename = 'TSNEPlot_aggr_filter_Spata22.tiff',width = 10, height =8)

p27<- FeaturePlot(aggr, features = "Ly6k", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p27
ggsave(filename = 'TSNEPlot_aggr_filter_Ly6k.tiff',width = 10, height =8)

p28<- FeaturePlot(aggr, features = "Tex12", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p28
ggsave(filename = 'TSNEPlot_aggr_filter_Tex12.tiff',width = 10, height =8)

p29<- FeaturePlot(aggr, features = "Ddx4", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p29
ggsave(filename = 'TSNEPlot_aggr_filter_Ddx4.tiff',width = 10, height =8)

p30<- FeaturePlot(aggr, features = "Gm4969", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p30
ggsave(filename = 'TSNEPlot_aggr_filter_Gm4969.tiff',width = 10, height =8)

p31<- FeaturePlot(aggr, features = "Tex12", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p31
ggsave(filename = 'TSNEPlot_aggr_filter_Tex12.tiff',width = 10, height =8)

p32<- FeaturePlot(aggr, features = "Ccna1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE)
p32
ggsave(filename = 'TSNEPlot_aggr_filter_Ccna1.tiff',width = 10, height =8)

p33<- FeaturePlot(aggr, features = "Ccnb3", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p33
ggsave(filename = 'TSNEPlot_aggr_filter_Ccnb3.tiff',width = 10, height =8)

p34<- FeaturePlot(aggr, features = "Rec8", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p34
ggsave(filename = 'TSNEPlot_aggr_filter_Rec8.tiff',width = 10, height =8)

p35<- FeaturePlot(aggr, features = "Sycp2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p35
ggsave(filename = 'TSNEPlot_aggr_filter_Sycp2.tiff',width = 10, height =8)

p36<- FeaturePlot(aggr, features = "Mei4", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.8)
p36
ggsave(filename = 'TSNEPlot_aggr_filter_Mei4.tiff',width = 10, height =8)

p37<- FeaturePlot(aggr, features = "Rnf212", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.8)
p37
ggsave(filename = 'TSNEPlot_aggr_filter_Rnf212.tiff',width = 10, height =8)

p38<- FeaturePlot(aggr, features = "Mnd1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.2)
p38
ggsave(filename = 'TSNEPlot_aggr_filter_Mnd1.tiff',width = 10, height =8)

p39<- FeaturePlot(aggr, features = "Mei1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.7)
p39
ggsave(filename = 'TSNEPlot_aggr_filter_Mei1.tiff',width = 10, height =8)

p40<- FeaturePlot(aggr, features = "Syce1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p40
ggsave(filename = 'TSNEPlot_aggr_filter_Syce1.tiff',width = 10, height =8)

p41<- FeaturePlot(aggr, features = "Syce2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2.2)
p41
ggsave(filename = 'TSNEPlot_aggr_filter_Syce2.tiff',width = 10, height =8)

p42<- FeaturePlot(aggr, features = "Rad51ap2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p42
ggsave(filename = 'TSNEPlot_aggr_filter_Rad51ap2.tiff',width = 10, height =8)

p43<- FeaturePlot(aggr, features = "Stag3", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.6)
p43
ggsave(filename = 'TSNEPlot_aggr_filter_Stag3.tiff',width = 10, height =8)

p44<- FeaturePlot(aggr, features = "4930447C04Rik", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p44
ggsave(filename = 'TSNEPlot_aggr_filter_4930447C04Rik.tiff',width = 10, height =8)

p45<- FeaturePlot(aggr, features = "Nipbl", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p45
ggsave(filename = 'TSNEPlot_aggr_filter_Nipbl.tiff',width = 10, height =8)

p46<- FeaturePlot(aggr, features = "Pttg1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.3)
p46
ggsave(filename = 'TSNEPlot_aggr_filter_Pttg1.tiff',width = 10, height =8)

p47<- FeaturePlot(aggr, features = "Rif1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.7)
p47
ggsave(filename = 'TSNEPlot_aggr_filter_Rif1.tiff',width = 10, height =8)

p48<- FeaturePlot(aggr, features = "Dnd1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p48
ggsave(filename = 'TSNEPlot_aggr_filter_Dnd1.tiff',width = 10, height =8)

p49<- FeaturePlot(aggr, features = "Mtor", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p49
ggsave(filename = 'TSNEPlot_aggr_filter_Mtor.tiff',width = 10, height =8)

p50<- FeaturePlot(aggr, features = "1700013H16Rik", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p50
ggsave(filename = 'TSNEPlot_aggr_filter_1700013H16Rik.tiff',width = 10, height =8)

p51<- FeaturePlot(aggr, features = "Taf7l", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p51
ggsave(filename = 'TSNEPlot_aggr_filter_Taf7l.tiff',width = 10, height =8)

p51<- FeaturePlot(aggr, features = "M1ap", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p51
ggsave(filename = 'TSNEPlot_aggr_filter_M1ap.tiff',width = 10, height =8)

p52<- FeaturePlot(aggr, features = "Hmgb2", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 3)
p52
ggsave(filename = 'TSNEPlot_aggr_filter_Hmgb2.tiff',width = 10, height =8)

p53<- FeaturePlot(aggr, features = "Smc3", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 2)
p53
ggsave(filename = 'TSNEPlot_aggr_filter_Smc3.tiff',width = 10, height =8)

p54<- FeaturePlot(aggr, features = "Parp1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p54
ggsave(filename = 'TSNEPlot_aggr_filter_Parp1.tiff',width = 10, height =8)

p55<- FeaturePlot(aggr, features = "Msh6", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p55
ggsave(filename = 'TSNEPlot_aggr_filter_Msh6.tiff',width = 10, height =8)

p56<- FeaturePlot(aggr, features = "Tex15", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1.5)
p56
ggsave(filename = 'TSNEPlot_aggr_filter_Tex15.tiff',width = 10, height =8)

p57<- FeaturePlot(aggr, features = "Tex14", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p57
ggsave(filename = 'TSNEPlot_aggr_filter_Tex14.tiff',width = 10, height =8)

p58<- FeaturePlot(aggr, features = "Ung", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p58
ggsave(filename = 'TSNEPlot_aggr_filter_Ung.tiff',width = 10, height =8)

p59<- FeaturePlot(aggr, features = "Samhd1", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p59
ggsave(filename = 'TSNEPlot_aggr_filter_Samhd1.tiff',width = 10, height =8)

p60<- FeaturePlot(aggr, features = "H1fx", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p60
ggsave(filename = 'TSNEPlot_aggr_filter_H1fx.tiff',width = 10, height =8)

p61<- FeaturePlot(aggr, features = "Msh5", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p61
ggsave(filename = 'TSNEPlot_aggr_filter_Msh5.tiff',width = 10, height =8)

p62<- FeaturePlot(aggr, features = "Ugt8a", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 0.5)
p62
ggsave(filename = 'TSNEPlot_aggr_filter_Ugt8a.tiff',width = 10, height =8)

p63<- FeaturePlot(aggr, features = "Xist", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, max.cutoff = 2)
p63
ggsave(filename = 'TSNEPlot_aggr_filter_Xist.tiff',width = 10, height =8)

p64<- FeaturePlot(aggr, features = "S100a4", reduction = "tsne", label = FALSE, pt.size = 2, cols = c("lightgrey", "red"), sort.cell = TRUE, min.cutoff = 1)
p64
ggsave(filename = 'TSNEPlot_aggr_filter_S100a4.tiff',width = 10, height =8)

aggr$seurat_clusters <- factor(aggr$seurat_clusters, levels = c('0', '2', '1', '3', '4'))
aggr@active.ident<- factor(aggr@active.ident, levels = c('0', '2', '1', '3', '4'))
new.cluster.ids <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4")
names(new.cluster.ids) <- levels(aggr)
aggr <- RenameIdents(aggr, new.cluster.ids)

colP<-c('#560047', '#a53bad', '#eb6bac', '#ffa8a0', '#f7dad7')
DimPlot(aggr, reduction = "tsne", label = FALSE, pt.size = 2, cols = colP)
ggsave(filename = 'TSNEPlot_aggr_filter_clusters.pdf',width = 10, height =8 )
colP2<-c('#560047', '#a53bad', '#eb6bac', '#ffa8a0', '#f7dad7')
DimPlot(aggr, reduction = "tsne", label = FALSE, pt.size = 2, group.by = 'TimePoints',cols = colP2)
ggsave(filename = 'TSNEPlot_aggr_filter_TimePoints.pdf',width = 10, height =8 )



VlnPlot(aggr, features = "Gfra1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Gfra1.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Lhx1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Lhx1.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Barhl2", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Barhl2.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Nefm", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Nefm.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Egr4", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Egr4.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Tcl1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Tcl1.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Nanos3", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Nanos3.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Etv5", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Etv5.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Bcl6b", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Bcl6b.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Zbtb16", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Zbtb16.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Stra8", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Stra8.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Rhox13", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Rhox13.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Kit", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Kit.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Sohlh1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Sohlh1.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Sohlh2", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Sohlh2.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Cyp26a1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Cyp26a1.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "Fbxo2", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Fbxo2.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Agpat3", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Agpat3.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Hist1h2aa", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Hist1h2aa.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Dmc1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Dmc1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Spo11", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Spo11.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Spata22", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Spata22.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Meiob", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Meiob.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Hormad1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Hormad1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Magea5", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Megea5.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Sycp3", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Sycp3.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Meioc", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Meioc.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Sycp1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Sycp1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Magea8", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Megea8.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Mybl1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Mybl1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Ly6k", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Ly6k.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Hist1h1a", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Hist1h1a.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Tex19.1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Tex19.1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Tslrn1", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Tslrn1.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Tex19.2", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Tex19.2.pdf',width = 10, height =8 )


VlnPlot(aggr, features = "Prdm9", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_Prdm9.pdf',width = 10, height =8 )

VlnPlot(aggr, features = "S100a4", slot = "counts", log = TRUE,cols = colP, pt.size = 0)
ggsave(filename = 'Violin plot_S100a4.tiff',width = 10, height =8 )


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
library(stringr)
s.genes=str_to_title(s.genes, locale = "en")
g2m.genes=str_to_title(g2m.genes, locale = "en")
aggr <- CellCycleScoring(aggr, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(aggr, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
aggr <- RunPCA(aggr, features = c(s.genes, g2m.genes))
aggr <- RunTSNE(aggr, dims = 1:30,features = c(s.genes, g2m.genes))
DimPlot(aggr, reduction = "tsne",group.by = "TimePoints")