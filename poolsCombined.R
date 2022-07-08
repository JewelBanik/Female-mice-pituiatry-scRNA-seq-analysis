setwd("C:/Users/JEWEL BANIK/Downloads/2021/UMAP_analysis_20genes")

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)

pool1.data <- Read10X(data.dir = paste(getwd(), "/pool1/", sep = ""))
pool1 <- CreateSeuratObject(counts = pool1.data, project = "pool1", min.cells = 3, min.features = 200)

pool2.data <- Read10X(data.dir = paste(getwd(), "/pool2/", sep = ""))
pool2 <- CreateSeuratObject(counts = pool2.data, project = "pool2", min.cells = 3, min.features = 200)

pool1
pool2

pools.combined <- merge(x = pool1, y=pool2, add.cell.ids=c("1", "2"), project="poolsCombined")

pools.combined
head(colnames(pools.combined))
table(pools.combined$orig.ident)
GetAssayData(pools.combined)[1:10, 1:10]


###QC
pools.combined[["percent.mt"]] <- PercentageFeatureSet(pools.combined, pattern = "^mt-")
head(pools.combined@meta.data, 5)

VlnPlot(pools.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pools.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pools.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

#### n.expressed.hkgenes <- ?(seurat@data[hkgenes.found, ] > ?)
#### seurat <- AddMetaData(object = ?, ? = ?, col.name = "n.exp.hkgenes")
n.expressed.hkgenes <- Matrix::colSums(pools.combined@assays$RNA['Gapdh', ] > 0)
pools.combined <- AddMetaData(object = pools.combined, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")
VlnPlot(object = pools.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","n.exp.hkgenes"), ncol = 4, pt.size = 0.1)
saveRDS(object = pools.combined, file = 'pools.combined.rawMerged.rds')

#discard unwanted cells that have >5% mitochondrial genes/counts and no. of Features <200 & >2500. 
# > summary(pools.combined@meta.data$nCount_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2794    5271    9116    9986   13032   94349 
# > IQR(pools.combined@meta.data$nCount_RNA)
# [1] 7760.5
# > 13032+(7760*1.5)
# [1] 24672
# > 5271-(7760*1.5)
# [1] -6369
# > summary(pools.combined@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05079 0.66601 0.90253 1.12795 1.33435 4.99347 
# > IQR(pools.combined@meta.data$percent.mt)
# [1] 0.6683385
# > 1.33435+(0.6683385*1.5)
# [1] 2.336858
# > 0.66601-(0.6683385*1.5)
# [1] -0.3364977
pools.combined <- subset(pools.combined.rawMerged, subset=nFeature_RNA>200&percent.mt<2.5&nCount_RNA<24700)

#Normalizing the data
pools.combined <- NormalizeData(pools.combined, normalization.method = "LogNormalize", scale.factor = 10000)

pools.combined <- FindVariableFeatures(pools.combined, selection.method = "vst", nfeatures = 2000)#vst: variance stabilizing transformation

###Scaling the Data: (all genes/features)
all.genes <- rownames(pools.combined)
pools.combined <- ScaleData(pools.combined, features = all.genes)
# pools.combined <- ScaleData(pools.combined, features = all.genes, vars.to.regress = 'percent.mt')

###Perform linear dimensional reduction: 
pools.combined <- RunPCA(pools.combined, features = VariableFeatures(object = pools.combined))
DimPlot(pools.combined, reduction = "pca")
VizDimLoadings(pools.combined, dims = 1:2, reduction = 'pca')
DimHeatmap(pools.combined, dims = 1, cells = 1000, balanced = T)
DimHeatmap(pools.combined, dims = 1:15, cells = 1000, balanced = T)

###Determine the 'dimensionality' of the dataset: 
pools.combined <- JackStraw(pools.combined, num.replicate = 100)
pools.combined <- ScoreJackStraw(pools.combined, dims = 1:20)
JackStrawPlot(pools.combined, dims = 1:20)

ElbowPlot(pools.combined)

saveRDS(object = pools.combined, file = 'pools.combined.jacked.rds')

###Cluster the cells: 
pools.combined15PCS<-FindNeighbors(pools.combined, dims=1:15)

pp <- FindClusters(object = pools.combined15PCS, resolution = 0.5)

pp <- RunUMAP(pp, dims = 1:15)
DimPlot(pp, reduction = 'umap', pt.size = 1, label = T)

table(pp@active.ident, pp@meta.data$orig.ident)
table(pp@meta.data$RNA_snn_res.0.5, pp@meta.data$orig.ident)
# https://github.com/satijalab/seurat/issues/738 
# https://github.com/satijalab/seurat/issues/2397

FeaturePlot(object = pp, features = c('Pou1f1', 'Gh', 'Ghrhr'), label = T)
FeaturePlot(object = pp, features = c('Pou1f1', 'Prl'), label = T, coord.fixed = T, min.cutoff = 'q85')
FeaturePlot(object = pp, features = c('Pou1f1', 'Tshb'), label = T)
FeaturePlot(object = pp, features = c('Fshb', 'Lhb', 'Gnrhr'), label = T)
FeaturePlot(object = pp, features = c('Pecam1'), label = T, coord.fixed = T)
FeaturePlot(object = pp, features = c('Pomc', 'Tbx19', 'Crhr1'), label = T, coord.fixed = T)
FeaturePlot(object = pp, features = c('Pou1f1', 'Prl', 'Gh', 'Top2a', 'Mki67'), label = F, coord.fixed = T)
FeaturePlot(object = pp, features = c('Hbb-bt'), label = T, coord.fixed = T, order = T)
FeaturePlot(object = pp, features = c('Pomc', 'Tbx19', 'Pax7'), label = T)
FeaturePlot(object = pp, features = c('Nkx2-1'), label = T, coord.fixed = T, order = T)
FeaturePlot(object = pp, features = c('C1qa'), label = T, order = T, coord.fixed = T)
FeaturePlot(object = pp, features = c('Col1a1'), label = T, order = T, coord.fixed = T)
FeaturePlot(object = pp, features = c('Sox2', 'Vim', 'Fstl1'), label = T, order = T, coord.fixed = T) #7.58%
FeaturePlot(object = pp, features = c('Msi1', 'Msi2', 'Prop1L'), label = T, coord.fixed = T, order = T)
FeaturePlot(object = pp, features = c('Atf4', 'Ormdl1', 'Ormdl3', 'Insig1'), label = T, 
            coord.fixed = T, min.cutoff = 0.85, order = T)
FeaturePlot(object = pp, features = c('Rab3d', 'Rab27a', 'Rab21', 'Rab2a'), label = T, 
            coord.fixed = T, min.cutoff = 0.85, order = T)



saveRDS(object = pp, file = "pp.UMAPed.rds", compress = 'gzip')

#resolution: 0.8, 1, and 1.2: 
pp2 <- FindClusters(object = pools.combined16PCS, resolution = c(0.8, 1, 1.2))

pp2 <- RunUMAP(pp2, dims = 1:16)
DimPlot(pp2, reduction = 'umap', pt.size = 1, label = T, group.by = 'RNA_snn_res.0.8')
DimPlot(pp2, reduction = 'umap', pt.size = 1, label = T, group.by = 'RNA_snn_res.1')
DimPlot(pp2, reduction = 'umap', pt.size = 1, label = T, group.by = 'RNA_snn_res.1.2')
DimPlot(pp, reduction = 'umap', pt.size = 1, label = T, group.by = 'RNA_snn_res.0.5')


FeaturePlot(object = pp2, features = c('Hbb-bt'), label = T, min.cutoff = 'q50')
FeaturePlot(object = pp2, features = c('Prl'), label = T, min.cutoff = 'q85')
VlnPlot(object = pp, features = c('Prl', 'Gh','Tshb', 'Pou1f1'), ncol = 1)
VlnPlot(object = pp, features = c('Hbb-bt', 'Prl'), ncol = 1)


#renaming clusters per cell types
# https://github.com/satijalab/seurat/issues/3239

# https://www.biostars.org/p/460935/#462724 #adding objects as metadata. 

FeaturePlot(pools.combined, features = c("Gh","Ghrhr", "Prl", "Tshb"))

FeaturePlot(pools.combined, features = c("Gnrhr","Lhb", "Fshb", "Nr5a1"))

FeaturePlot(pools.combined, features = c("Pomc","Pax7", "Tbx19"))

FeaturePlot(pools.combined, features = c("Pou1f1","Sox2", "Sox9","Mki67"))
FeaturePlot(pools.combined, features = c("Prop1L","Nkx2-1", "Pecam1","Col1a1", 'Hbb-bt', 'C1qa'))


FeaturePlot(pools.combined, features = c("Norad", "Pitx1", "Lsm14a","Lsm14b" ,"Cldn4", "Cxcl12"))

FeaturePlot(pools.combined, features = c("Foxp2", "Aif1l", "Ihh","Avpr1b" ,"Ghsr", "Foxo3"))

FeaturePlot(pools.combined, features = c("Atf4", "Rxra", "Ddx6","Foxj1"))

FeaturePlot(pools.combined, features = c("Nkx2-1", "Pecam1", "Col1a1","Hbb-bt" ,"C1qa"))

FeaturePlot(pools.combined, features = c("Msi1", "Msi2", "Prop1L","Sox2"))

###Assigning cell type identity to clusters: 
new.cluster.ids <- c("Lactotropes", "Lactotropes", "Lactotropes" , "Lactotropes", "Somatotropes", "Lactotropes", "Pituitary Stem Cells","Corticotropes","Melanotropes", "Pituitary Stem Cells", "Endothelia", "Gonadotropes", "WBC", "Posterior Pituitary", "Thyrotropes", "Connective Tissue", "Proliferating Cells")
names(new.cluster.ids) <- levels(pp)
pit.clustersNames <- RenameIdents(pp, new.cluster.ids)
DimPlot(pit.clustersNames, reduction = "umap", label = TRUE, pt.size = 0.5)+coord_fixed()
pit.clustersNames$cell_types <- pit.clustersNames@active.ident
View(pit.clustersNames@meta.data)
saveRDS(object = pit.clustersNames, file = 'pit.ClustersNames.final.rds', compress = T)

idOrders <- c("Somatotropes", "Lactotropes","Thyrotropes", "Gonadotropes", "Corticotropes","Melanotropes", "Pituitary Stem Cells", "Proliferating Cells", "Posterior Pituitary", "Endothelial Cells", "WBC", "Connective Tissue")
# Rename identity classes
#pit.clustersNames <- RenameIdents(object = pit.clustersNames, `Endothelial Cells` = "Endothelia")
#pit.clustersNames$cell_types <- pit.clustersNames@active.ident


###Dim and Feature plots: 
png(file = "UMAP_cell_types.png", height = 4, width = 8, res = 1200, units = 'in')
DimPlot(pit.clustersNames, reduction = 'umap', group.by = 'cell_types', label = T,
        repel = F, label.size = 2, pt.size = 0.25)+coord_fixed()+NoLegend()
dev.off()

png(file = "FeaturePlotCombined.png", height = 7, width = 12, res = 1200, units = 'in')
plot1 <- DimPlot(object = pit.clustersNames, reduction = 'umap', label = T, pt.size = 1,
                 repel=F, label.size=2, group.by = 'cell_types', combine = T)& 
  coord_equal(ratio = 1)& labs(title = 'Female Mice Pituitary')&NoLegend()
plot2 <- FeaturePlot(object = pit.clustersNames, features = 'Ormdl1', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot3 <- FeaturePlot(object = pit.clustersNames, features = 'Ormdl3', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot4 <- FeaturePlot(object = pit.clustersNames, features = 'Insig1', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

png(file = "FeaturePlotCombined2.png", height = 7, width = 12, res = 1200, units = 'in')
plot1 <- DimPlot(object = pit.clustersNames, reduction = 'umap', label = T, pt.size = 1,
                 repel=T, label.size=2, group.by = 'cell_types', combine = T)& 
  coord_equal(ratio = 1)& labs(title = 'Female Mice Pituitary')&NoLegend()
plot2 <- FeaturePlot(object = pit.clustersNames, features = 'Rab2a', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot3 <- FeaturePlot(object = pit.clustersNames, features = 'Rab3d', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot4 <- FeaturePlot(object = pit.clustersNames, features = 'Rab21', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot5 <- FeaturePlot(object = pit.clustersNames, features = 'Rab27a', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)

CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), ncol = 3)

dev.off()


png(file = "FeaturePlotCombined3.png", height = 7, width = 10, res = 1200, units = 'in')
plot1 <- DimPlot(object = pit.clustersNames, reduction = 'umap', label = T, pt.size = 1,
                 repel=T, label.size=2, group.by = 'cell_types', combine = T)& 
  coord_equal(ratio = 1)& labs(title = 'Female Mice Pituitary')&NoLegend()
plot6 <- FeaturePlot(object = pit.clustersNames, features = 'Atf4', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot7 <- FeaturePlot(object = pit.clustersNames, features = 'Creb3l2', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
CombinePlots(plots = list(plot1, plot6, plot7), ncol = 2)

dev.off()

png(file = "FeaturePlotCombined4.png", height = 7, width = 10, res = 1200, units = 'in')
plot1 <- DimPlot(object = pit.clustersNames, reduction = 'umap', label = T, pt.size = 1,
                 repel=T, label.size=2, group.by = 'cell_types', combine = T)& 
  coord_equal(ratio = 1)& labs(title = 'Female Mice Pituitary')&NoLegend()
plot6 <- FeaturePlot(object = pit.clustersNames, features = 'Msi1', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot7 <- FeaturePlot(object = pit.clustersNames, features = 'Msi2', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
plot8 <- FeaturePlot(object = pit.clustersNames, features = 'Prop1L', cols = c('grey', 'blue'), 
                     order = T, coord.fixed = T, combine = T, min.cutoff = 'q80', pt.size = 1)
CombinePlots(plots = list(plot1, plot6, plot7, plot8), ncol = 2)

dev.off()


###Dot Plot: 
png(file = "DotPlot_5pctFiltered.png", height = 5, width = 10, res = 1200, units = 'in')
DotPlot(pit.clustersNames, features = c("Msi1","Msi2", "Ormdl1", 'Ormdl3', 'Insig1', 'Rab2a', 
                           'Rab3d', 'Rab21', 'Rab27a', 'Atf4', 'Creb3l2'), cols = c('grey', 'red'),
        group.by = "cell_types", dot.min = 0.05)+labs(title = "Mouse Pituitary")+theme_bw()
dev.off()

###Violin plots: 
features <- c("Msi1","Msi2", "Ormdl1", 'Ormdl3', 'Insig1', 'Rab2a', 
              'Rab3d', 'Rab21', 'Rab27a', 'Atf4', 'Creb3l2')
png(file = "VlnPlot_cell_types.png", height = 7, width = 10, res = 1200, units = 'in')
violin1 <- VlnPlot(object = pit.clustersNames, features = features, group.by = 'cell_types', 
                   ncol = 1, pt.size =1.5, combine = T, same.y.lims = T, stack = T, 
                   flip = T)& 
  NoLegend()& labs(title = "Mouse Pituitary")

violin1
dev.off()


# to find out how many cells in a cluster express a gene: 
Rab27a <- table(pit.clustersNames@assays$RNA@data['Rab27a', ]>0, pit.clustersNames@meta.data$cell_types)
Rab27a <- as.data.frame(t(Rab27a))

write.table(x = Rab27a, file = "Rab27a.txt", quote = T, sep = "\t", row.names = T, col.names = T)

library(stargazer)
stargazer(as.data.frame(t(Rab27a)), type = "text", summary = F, out = "Rab27a_PerCluster.txt", 
          colnames = T)

lac <- function(y){
  y = "Lactotropes"
  lactotropes <-  Rab27a$Freq[grep(pattern = y , x = Rab27a$Var1, ignore.case = T, value = F)]
  lac <- (lacto[2]/sum(lacto))*100
  paste(lac,'%', sep = "")
}


Rab27a$Freq[grep(pattern = 'Lactotropes', x = Rab27a$Var1, ignore.case = T, value = F)]
lacto <- Rab27a[ , ]$Freq[c(1, 13)]
lacto <- (lacto[2]/sum(lacto))*100
paste(lacto,'%', sep = "")


# i= 0
# for (i in seq(1,10)){
#   if (i==4)
#     next
#  print(i) 
# }

# Plot order VlnPlot #454
# How can I change the x-axis ordering in VlnPlot? #1250
###to have violin plot with sorted cell types: load the .rds file first. 
# Define an order of cluster identities
levels(pit.clustersNames)
my_levels <- c('Somatotropes', 'Lactotropes', 'Thyrotropes', 'Corticotropes', 
               'Gonadotropes', 'Pituitary Stem Cells', 'Proliferating Cells',
               'Melanotropes', 'Posterior Pituitary', 'Endothelia', 
               'Connective Tissue', 'WBC')

# Relevel object@ident
pit.clustersNames@active.ident <- factor(x = pit.clustersNames@active.ident,
                                               levels = my_levels)


features <- c("Msi1","Msi2", "Ormdl1", 'Ormdl3', 'Insig1', 'Rab2a', 
              'Rab3d', 'Rab21', 'Rab27a', 'Atf4', 'Creb3l2')
png(file = "VlnPlot_cell_types_sorted.png", height = 7, width = 10, res = 1200, units = 'in')
violin1 <- VlnPlot(object = pit.clustersNames, features = features, pt.size =1.5, 
                   combine = T, same.y.lims = T, stack = T, flip = T)& 
  NoLegend()& labs(title = "Female Mouse Pituitary")

violin1
dev.off()


