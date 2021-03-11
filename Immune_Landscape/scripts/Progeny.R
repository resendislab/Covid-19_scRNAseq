rm(list = ls())

library(progeny)
library(Seurat)
library(pheatmap)


## i is replaced for any of the following cell types:
######## Monocytes
######## Macrophages
######## Dendritic
######## T.NK  
######## B
######## Ephitelial

####### For instance
i <- "Dendritic"
### Loading data
balf<-readRDS(paste("./Results/Data/Clusters/balf.",
                    i,".rds",sep = ""))

########### Seurat standard workflow
balf <- RunPCA(balf, npcs = 30, verbose = TRUE)
### t-SNE and Clustering
set.seed(0)
balf <- RunUMAP(balf, reduction = "pca", dims = 1:20)
# set.seed(0)
balf <- FindNeighbors(balf, reduction = "pca",
                      dims = 1:20)
balf<- FindClusters(balf, resolution = 0.5)

DimPlot(balf, reduction = "umap")

#### Number of clusters
n <- max(as.integer(as.vector(balf@meta.data$seurat_clusters)))
new.cluster.ids <- as.character(0:n)
names(new.cluster.ids) <- levels(balf)
balf <- RenameIdents(balf, new.cluster.ids)

## Data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
CellsClusters <- data.frame(Cell = names(Idents(balf)), 
                            CellType = as.character(Idents(balf)),
                            stringsAsFactors = FALSE)
###Progeny computed 
## Progeny activity scores and their addition to the Seurat object
## as a new assay called Progeny. 
balf <- progeny(balf, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## Seurat functions applied to the Progeny scores. 
## For instance, pathway activity scores scalation. 
balf <- Seurat::ScaleData(balf, assay = "progeny") 

## Progeny scores trnasformation into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(balf, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## Matching Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

##### Add column organized by health status (Svr|Crt|Mld)
progeny_scores_df$Health_status <- gsub("(Svr|Crt|Mld).*", "\\1", 
                                        progeny_scores_df$Cell)
##### Add column organized by health status (Svr|Crt|Mld) and cluster assignation
progeny_scores_df$Mix <- paste(progeny_scores_df$Health_status,
                               progeny_scores_df$CellType, sep = "_")
#####

progeny_scores_df$Health_status <- gsub("Svr","Severe",progeny_scores_df$Health_status)
progeny_scores_df$Health_status <- gsub("Mld","Moderate",progeny_scores_df$Health_status)
progeny_scores_df$Health_status <- gsub("Crt","Control",progeny_scores_df$Health_status)

nms <-  c("cell_type","Health_status","Mix")
for(j in 1:3){
  if(j == 1){
    ## Progeny scores by cellpopulation
    summarized_progeny_scores <- progeny_scores_df %>% 
      group_by(Pathway, CellType) %>%
      summarise(avg = mean(Activity), std = sd(Activity))
  }else if(j == 2){
    ## Progeny scores by Health_Status
    summarized_progeny_scores <- progeny_scores_df %>% 
      group_by(Pathway, Health_status) %>%
      summarise(avg = mean(Activity), std = sd(Activity))
  }else {
    ## Progeny scores by clustering and health status
    summarized_progeny_scores <- progeny_scores_df %>% 
      group_by(Pathway, Mix) %>%
      summarise(avg = mean(Activity), std = sd(Activity))
  }
  
  ###Plot progeny 
  # Data preparation for the plot
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

  breaksList = seq(-1.1, 1.1, by = .05)
  myColor =colorRampPalette(c("Darkblue", "white","red"))(length(breaksList))
  progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),
                         breaks = breaksList,
                         color = myColor,
                         main = "PROGENy", 
                         angle_col = 0, 
                         treeheight_col = 2,  border_color = NA)
}
