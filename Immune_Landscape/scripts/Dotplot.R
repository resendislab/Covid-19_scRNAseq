rm(list = ls())

library(Seurat)

## i is replaced for any of the following cell types:
######## Monocytes
######## Macrophages
######## Dendritic
######## T.NK  
######## B
######## Ephitelial

i <- "Monocytes"


#### Loading Data
balf<-readRDS(paste("./Results/Data/Clusters/pbmc.",
                    i,".rds",sep = ""))

genes = c("CD14","FCGR3A","MRC1","MME","CXCR3","IL3RA","CD1C",
          "HLA-DRB1","HLA-DRA","IFI6","ISG15","CXCL11","CCR2",
          "CD68","CD69","CXCL10")
balf$orgn <- gsub("Mld","Mod", balf$orgn) 

features <- genes[(which(genes %in% rownames(balf@assays$RNA)))]
balf.Dot <- DotPlot(balf, features = genes, 
                    cols = "RdYlBu",group.by = "orgn",
                    assay = "RNA",dot.scale = 4)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  RotatedAxis()+
  theme(text = element_text(size = 11),
        axis.text = element_text(size = 9))+
  guides(colour = guide_colourbar(title = 'Average\nexpression'),
         size=guide_legend(title = "% Expressed"))
print(balf.Dot)

