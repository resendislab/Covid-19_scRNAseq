#### This scripts downloads data directly from GEO, the used repository are
#### related to the paper paper https://www.nature.com/articles/s41591-020-0901-9.

rm(list = ls())
setwd("./Covid-19_scRNAseq")

library(GEOquery)
library(stringr)

#### Function to download data
Download <- function(GEO.acs.nbr){
  nms <- names(GEO.acs.nbr)
  for (i in nms){ 
    pth <- paste("./Data/",i,sep = "")
    if (dir.exists(pth)==FALSE){dir.create(pth)}
    lapply(GEO.acs.nbr[[1]][1], 
                     function(x) getGEOSuppFiles(x,baseDir = pth, 
                                                 makeDirectory = FALSE))
  }
}
#### Folder to save data
if (dir.exists("./Data")==FALSE){dir.create("./Data")}
#### Accesion information
fl <- getGEO(GEO = "GSE145926", destdir = "./Data")
#### samples number
GEO.acs <- fl$GSE145926_series_matrix.txt.gz$geo_accession

#### Severe
Svr <- GEO.acs[str_detect(GEO.acs,"(69|70|71|72|73|74)$")]
#### Mild
Mld <- GEO.acs[str_detect(GEO.acs,"(48|49|50)$")]
#### healthy Control
Crt <- GEO.acs[str_detect(GEO.acs,"(51|52|53)$")]
GEO.acs.nbr <- list("Severe"=Svr,"Mild"=Mld, "Control"=Crt)

Download(GEO.acs.nbr)

