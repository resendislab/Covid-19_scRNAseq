
rm(ls=list())
library(Seurat)

### Loading data processed on the Seurat.R file
balf.combined <- readRDS(file = "./Data/Results/balf.rds")

### Association between cell types and clusters
clust.idents <- list(
  Monocytes = c(10),
  Macrophages = c(0:4,6,11,19),
  Dendritic = 7,
  T.NK = c(8,9,13,14,16),  
  B = c(12,18),
  Ephitelial = c(15,17)
  #  Other = 5 
)

##### By Identity
nm <- names(clust.idents)
invisible(lapply(1:length(clust.idents),
                 function(x){
                   ### Separating celltypes data
                   balf1 <- subset(balf.combined, idents = as.character(x))
                   ### Saving cell types data
                   saveRDS(balf1, paste("./Results/Data/Clusters/balf.",
                                        nm[x],".rds", sep = ""))
                   }))