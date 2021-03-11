library(RColorBrewer)
library(ggplot2)

umap_info <- read.csv("./Results/Data/Clusters/UMAP_info.csv")

clust.idents <- list(
  Monocytes = c(10),
  Macrophages = c(0:4,6,11,19),
  Dendritic = 7,
  T.NK = c(8,9,13,14,16),  
  B = c(12,18),
  Ephitelial = c(15,17),
  Other = 5 
)

Pheno <- vector("character", length = length(umap_info$clust_info))
Pheno.colors<- brewer.pal(n = length(clust.idents), name = "Set1")
Prop <- matrix(0,3,length(clust.idents), 
                      dimnames = list(c("Control","Moderate","Severe"),names(clust.idents)))

for(x in 1:length(clust.idents)){
  ### Association between samples and cell type
  idx = which(umap_info$clust_info %in% clust.idents[[x]])
  Pheno[idx] = names(clust.idents)[x]
  #### Proportion for each cell type scaled for the number of patients for each group
  Prop[,x] <- table(gsub("(Svr|Mld|Crt).*","\\1",rownames(umap_info)[idx]))/c(3,3,6)
  #### Scaling for the number of samples in every celltype
  Prop[,x] <- Prop[,x]/sum(Prop[,x])*100
}

## uMAP visualization
umap_info <- cbind(umap_info, Pheno)
umap_info$Pheno <- gsub("T.NK","T & NK",umap_info$Pheno)
nms <- c("Monocytes","Macrophages","Dendritic","T & NK","B","Ephitelial","Other")

umap_tot <- ggplot(umap_info, aes(x=UMAP_1, y=UMAP_2, fill= factor(Pheno, level = nms)))+ 
  labs(x= "uMAP 1",y ="uMAP 2")+
  geom_point(size=2, shape=21)+ 
  theme(legend.title = element_blank(),legend.background=element_blank())+
  scale_fill_manual(values = Pheno.colors)+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  theme(text = element_text(size = 11),axis.text = element_text(size = 9))+
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())+
  theme(legend.key.height=unit(0.07,"line"),
        legend.justification ="top",legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-10,-5))
print(umap_tot)

## barplot visualization
ZZ <- data.frame(HS=rownames(Prop)[row(Prop)], 
                 Cell_type=colnames(Prop)[col(Prop)],
                 Proportion=c(Prop))

ZZ$HS <- factor(ZZ$HS , levels = c("Severe","Moderate","Control"))
ZZ$Cell_type <- gsub("T.NK", "T & NK",ZZ$Cell_type)

barplot_prop <- ggplot(ZZ, aes(x= Proportion,fill=HS, y=  factor(Cell_type, level = rev(unique(Cell_type))))) + 
  geom_bar(stat = "identity",position=position_dodge())+
  scale_fill_manual(values=rev(c(wes_palette(n=3, name="Darjeeling1"))),
                    guide = guide_legend(reverse=TRUE),
                    labels = c( "Severe", "Moderate","Control"))+
  theme_bw()+xlab("Proportion (%)") +
  theme(legend.position="bottom",
        axis.title.y = element_blank(),
        legend.title=element_blank())+
  theme(text = element_text(size = 11),
        axis.text = element_text(size = 9))
print(barplot_prop)
