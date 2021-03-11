rm(list = ls())

##### 
## Required packages
library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

## i is replaced for any of the following cell types:
######## Monocytes
######## Macrophages
######## Dendritic
######## T.NK  
######## B
######## Ephitelial
## For instance
i <- "Dendritic"

#### Loading Data
balf<-readRDS(paste("./Results/Data/Clusters/balf.",
                    i,".rds",sep = ""))
#########################################
#########################################

########### Seurat standard workflow
balf <- RunPCA(balf, npcs = 30, verbose = TRUE)
### t-SNE and Clustering
set.seed(0)
balf <- RunUMAP(balf, reduction = "pca", dims = 1:20)
# set.seed(0)
balf <- FindNeighbors(balf, reduction = "pca",
                       dims = 1:20)
balf<- FindClusters(balf, resolution = 0.5)

##### Number of cluster extracted from the seurat object
n <- max(as.integer(as.vector(balf@meta.data$seurat_clusters)))

## Assigning cell type identity to clusters
new.cluster.ids <- as.character(0:n)
names(new.cluster.ids) <- levels(balf)

balf <- RenameIdents(balf, new.cluster.ids)

DimPlot(balf, reduction = "umap", label = TRUE, pt.size = 0.5)+ NoLegend()

## Regulons based on interactions with confidence level A, B and C
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

balf <- run_viper(balf, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))

## `r CRANpkg("Seurat")` to cluster the cells following the same protocol than above but using TF activity scores.

## Nearest Neighbours computation to perform cluster based on TF
DefaultAssay(object = balf) <- "dorothea"
balf <- ScaleData(balf)
balf <- RunPCA(balf, features = rownames(balf),
                verbose = FALSE)
balf <- FindNeighbors(balf, dims = 1:50, verbose = FALSE)
balf <- FindClusters(balf, resolution = 0.8, 
                      verbose = FALSE)
balf <- RunUMAP(balf, dims = 1:50, 
                 umap.method = "uwot", metric = "cosine")
#balf.markers.TF <- FindAllMarkers(balf, only.pos = TRUE, 
#                                   min.pct = 0.25, 
#                                   logfc.threshold = 0.25, 
#                                   verbose = FALSE)

##### Number of cluster
n <- max(as.integer(as.vector(balf@meta.data$seurat_clusters)))
## Assigning cell type identity to clusters
new.cluster.ids <- as.character(0:n)

names(new.cluster.ids) <- levels(balf)
balf <- RenameIdents(balf, new.cluster.ids)

DimPlot(balf, reduction = "umap", 
        label = TRUE, pt.size = 0.5)

## TF activity per cell population
## Viper scores transformation, scaled by seurat, into a data frame to better 
## handling the results

viper_scores_df <- GetAssayData(balf, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame() %>%
  t()
#write.csv(viper_scores_df, "./Results/Dorothea/ViperScoresDendritic.csv")

## Data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(balf)), 
                            cell_type = as.character(Idents(balf)),
                            stringsAsFactors = FALSE)

## Data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

##### Adding column by health status (Svr|Crt|Mld)
viper_scores_clusters$Health_status <- gsub("(Svr|Crt|Mld).*", "\\1", 
                                            viper_scores_clusters$cell)

viper_scores_clusters$Health_status <- gsub("Svr", "Severe", 
                                            viper_scores_clusters$Health_status)
viper_scores_clusters$Health_status <- gsub("Mld", "Moderate", 
                                            viper_scores_clusters$Health_status)
viper_scores_clusters$Health_status <- gsub("Crt", "Crontol", 
                                            viper_scores_clusters$Health_status)

##### Adding column by health status (Svr|Crt|Mld) and cluster assignation
viper_scores_clusters$Mix <- paste(viper_scores_clusters$Health_status,
                                   viper_scores_clusters$cell_type, sep = "_")
#####

#### Cluster
nms <-  c("cell_type","Health_status","Mix")
for(j in 1:3){
  # Cell population
  if(j == 1){
    summarized_viper_scores <- viper_scores_clusters %>%
      group_by(tf, cell_type) %>%
      summarise(avg = mean(activity),
                std = sd(activity))
  }else if(j ==2){
    #### Health Status
    summarized_viper_scores <- viper_scores_clusters %>% 
      group_by(tf, Health_status) %>%
      summarise(avg = mean(activity),
                std = sd(activity))
  }else{
    #### Mix
    summarized_viper_scores <- viper_scores_clusters %>% 
      group_by(tf, Mix) %>%
      summarise(avg = mean(activity),
                std = sd(activity))
  }
  ## Selecction of the 50 most variable TFs.
  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(50*n, var) %>%
    distinct(tf)

  ## data preparation for the plot
  summarized_viper_scores_df <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%   
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, 
               stringsAsFactors = FALSE) 

  breaksList = seq(-1.1, 1.1, by = .05)
  my_color=colorRampPalette(c("Darkblue", "white","red"))(length(breaksList))
  viper_hmap = pheatmap(t(summarized_viper_scores_df),
                         color=colorRampPalette(c("Darkblue", "white","red"))(length(breaksList)),
                         breaks = breaksList,
                         main = "DoRothEA",
                         angle_col = 0,
                         treeheight_col = 2,
                         border_color = NA,fontsize_row = 9)
}
