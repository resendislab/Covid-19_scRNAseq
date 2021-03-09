### 

rm(list = ls())
setwd("./Covid-19_scRNAseq")

load.sparse <- function(pct,nms){
  fls <- list.files(paste("./Data/",pct,sep = ""),"*.h5$")
  
  h5.fls <- lapply(fls, function(x) read.and.labels(pct, x, nms))
  
  #### Merged files
  h5.merged <- merge.sparse(h5.fls)
  
  return(h5.merged)
}

read.and.labels = function(pct,x,nms) {
  data <- Read10X_h5(paste("./Data",pct,x,sep = "/"))
  lbl <- str_extract(x,"C\\d*")
  nms.fl <- unlist(lapply(1:length(data@Dimnames[[2]]),
                           function(x) paste(nms,lbl,x,sep = "_")))
  data@Dimnames[[2]] = nms.fl
  
  return(data)
}

merge.sparse = function(M.list) {
  A = M.list[[1]]
  if (length(M.list) > 1){
    for (i in 2:length(M.list)){ #i indexes of matrices
      # finding what's missing
      misA = rownames(M.list[[i]])[!rownames(M.list[[i]]) %in% 
                                     rownames(A)]
      misB = rownames(A)[!rownames(A) %in% rownames(M.list[[i]])]
      
      misAl = as.vector(numeric(length(misA)), "list")
      names(misAl) = misA
      misBl = as.vector(numeric(length(misB)), "list")
      names(misBl) = misB
      
      ## adding missing columns to initial matrices
      An = Reduce(rbind, c(A, misAl))
      if (length(misA) > 0)
      {
        lenA <- nrow(An)-length(misA)+1
        rownames(An)[lenA:nrow(An)] = names(misAl)
      }
      
      Bn = Reduce(rbind, c(M.list[[i]], misBl))
      if(length(misB) > 0)
      {
        lenB <- nrow(Bn)-length(misB)+1
        rownames(Bn)[lenB:nrow(Bn)] = names(misBl)
      }
      
      Bn <- Bn[rownames(An),]
      
      # final bind
      A = cbind(An, Bn)
      #    print(c(length(M.list), i))
    }
  }
  return(A)
}

library(rhdf5)
library(Seurat)
library(Matrix)
library(stringr)

####### Raw sparce Matrix dgCMatrix
## Severe data files
Svr.merged <- load.sparse("Severe","Svr")
## Moderate data files
Mld.merged <- load.sparse("Mild","Mld")
## Moderate data files
Crt.merged <- load.sparse("Control","Crt")
##### Data with all experimental groups
MyData <- merge.sparse(list(Svr.merged,Mld.merged,Crt.merged))

rm(Crt.merged,Mld.merged,Svr.merged)

if (dir.exists("./Data/Processed")==FALSE){dir.create("./Data/Processed")}
writeMM(MyData,"./Data/Processed/Data.txt")

### Filtered Matrix dgCMatrix, samples with less than nonzeros values
### <=200 are elimanited and genes into less than 3 samples are elimanted  

metaData<- data.frame(
  #### origin (Severe, mild and Control)
  orgn <- gsub("[:_:]C\\d*[:_:][0-9]*$","",colnames(MyData)),
  #### sample
  smpl <- str_extract(colnames(MyData),"C\\d{2,3}")
)
rownames(metaData) <- colnames(MyData)
colnames(metaData) <- c("orgn","smpl")

MyData.1 <- CreateSeuratObject(counts = MyData, 
                               min.cells = 3, 
                               min.features = 200, 
                               project = "test",
                               names.field = 1,
                               names.delim = "_",
                               meta.data = metaData)

rm(MyData)

###### Saving Data for the immune exploration analysis
saveRDS(object = MyData.1, file = "./Data/MyDataCovid.rds")

MyData.2 <- MyData.1@assays$RNA@counts

###### Saving Data for the machine learning analysis
## Count Matrix
writeMM(MyData.2,"./Data/Processed/Data_filtered.txt")
## Genes
write(x = rownames(MyData.2), file = "./Data/Processed/genes.tsv")
## Samples
write(x = colnames(MyData.2), file = "./Data/Processed/labels.tsv")
