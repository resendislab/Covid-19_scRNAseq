### This script perform seurat implementation: Data were analyzed to get filtering, clustering and Diffenriatial Gene Expresion analyses.

rm(list = ls())

setwd("./Covid-19_scRNAseq")

library(Seurat)
library(ggplot2)
library(tibble)
library(cowplot)

MyData.1 <- readRDS(file = "./Data/MyDataCovid.rds")


## Identification of mitochondrial genes
MyData.1[["percent.mt"]] <- PercentageFeatureSet(MyData.1, pattern = "^MT-")
MyData.1 <- subset(MyData.1, percent.mt < 5)

balf <- SplitObject(MyData.1, split.by = "orgn")
rm(MyData.1)

####### Pre-processing, normalization and identification of highly variable features

balf.list <- lapply(X = balf, FUN = function(x) {
 x <- NormalizeData(x, normalization.method = "LogNormalize",scale.factor = 10000)
 x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

balf.anchors <- FindIntegrationAnchors( object.list = balf.list)
rm(balf.list,balf)
balf.combined <- IntegrateData(anchorset = balf.anchors)
rm(balf.anchors)
DefaultAssay(balf.combined) <- "integrated"

### Run the standard workflow for visualization and clustering
set.seed(0)
balf.combined <- ScaleData(balf.combined, verbose = TRUE)
balf.combined <- RunPCA(balf.combined, npcs = 30, verbose = TRUE)
### uMAP and Clustering
balf.combined <- RunUMAP(balf.combined, reduction = "pca", dims = 1:20)
balf.combined <- FindNeighbors(balf.combined, reduction = "pca", 
                               dims = 1:20)
balf.combined <- FindClusters(balf.combined, resolution = 0.5)
balf.combined  <- RunUMAP(balf.combined, dims = 1:10,  umap.method = 'uwot', metric='cosine')

################# Saving umap projection coordinates and cluster information
if (dir.exists("./Results/Data/Clusters")==FALSE){dir.create("./Results/Data/Clusters",recursive = TRUE)}
Umap.info <- cbind(pbmc.combined@reductions$umap@cell.embeddings,
                   as.integer(as.vector(pbmc.combined@meta.data$seurat_clusters)))

write.csv(Umap.info,"./Results/Data/Clusters/UMAP_info.csv")

################### visualization
if (dir.exists("./Results/Plots")==FALSE){dir.create("./Results/Plots",recursive = TRUE)}

png(file='./Results/Plots/umap.png', width =700, 
    height = 700, units = "px")
p1 <- DimPlot(pbmc.combined, reduction = "umap", label = TRUE,
              label.size = 8) + NoLegend()+ labs(x="uMAP 1", y= "uMAP 2")
plot_grid(p1)
dev.off()
#################

## Finding differentially expressed features (cluster biomarkers)
balf.markers <- FindAllMarkers(balf.combined, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, verbose = FALSE)
pth <- "./Results/Data/Clusters"
if (dir.exists(pth)==FALSE){dir.create(pth,recursive = TRUE)}

invisible(lapply(0:n, function(x) write.table(
  balf.markers$gene[which(balf.markers$cluster == x)],
  file = paste(pth,"/DGA_cluster_ALL_",x,".txt", sep = ""),
  quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "/t")))

saveRDS(object = balf.combined, file = "./Data/balf.rds")