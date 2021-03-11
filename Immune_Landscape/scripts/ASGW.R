
rm(list = ls())

library(readxl)
library(purrr)
######
list.as.matrix <- function(x, byrow=FALSE, filler=""){
  maxlen <- max(sapply(x,length))
  xm <- sapply(x, 
               function(xs){
                 if (!identical(unlist(MakGns1), character(0))){
                   fillen <- maxlen-length(xs)
                   if (fillen > 0) {
                     c(xs,rep(filler,fillen))
                   } else xs
                 } else filler })
  if (byrow)return(t(xm)) else return(xm)
}
###### Read LM22 marker genes xls file manually curated
sheets <- excel_sheets("./Data/ASGW.xlsx")
ASGW <- map_df(sheets, ~ read_excel("./Data/ASGW.xlsx", sheet = .x, col_types = "text"))

n <- dim(ASGW)[1]
lbs <- vector("character",length = n)
#blk.spc <- vector("integer",0)
spc.gns <- vector("integer",0)
for (i in 1:n){
  flg <- ASGW$`Cell Type`[i]
  if (!is.na(flg)){ flg2 = flg }
  lbs[i] <- flg2
  if ((ASGW$ASGW[i] == "Specific" || 
       ASGW$ASGW[i] == 100 || ASGW$ASGW[i] == 75 ||
      ASGW$ASGW[i] == "Partially specific") &&
      !is.na(ASGW$ASGW[i])){
    spc.gns <- c(spc.gns,i)
  }
}
ASGW$`Cell Type` <- lbs
#### Genes with the label "Specific" and "Partially Specific"
#### Blank rows were elimited and Cell Type information was completed
ASGW1 <- ASGW[spc.gns,]

ASGW2 <- split(ASGW1$Gene, ASGW1$`Cell Type`)

### Load marker genes for every cluster
mks <- lapply(0:19, function(x) 
  read.table(
    paste("./Results/Data/Clusters/DGA_cluster_ALL_",x,
          ".txt", sep = ""), header = FALSE, col.names = FALSE, 
    sep = "\t"))

names(mks) <- unlist(lapply(0:(length(mks)-1), function(x) 
  paste("Cluster",x, sep = "_")))

b <- matrix(data = NA, nrow = 1, ncol = 21)
colnames(b) <- c("V1", names(mks))


### Comparing both list
for (i in names(ASGW2)){
  MakGns <- ASGW2[[i]]
  MakGns1 <- lapply(mks, function(x) 
    MakGns[MakGns %in% unlist(x)])
  a <- list.as.matrix(MakGns1)
#  print(a)
  n <- length(a)
  if (n > 20) {  
    a <- cbind(rep(i,nrow(a)),a)
  }else{
      a <- c(i,a)
  }
  b <- rbind(b,a)
}

b <- b[-1,]

if (dir.exists("./Results/Data/ASGW")==FALSE){
  dir.create("./Results/Data/ASGW",recursive = TRUE)}
write.table(b,"./Results/Data/ASGW/Summary.csv", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")





