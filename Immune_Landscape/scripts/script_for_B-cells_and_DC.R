## WebGestaltR
# load library
library(WebGestaltR)
# Working directory
Folder <- "~/pathway/to/working/directory/"
# Output directory
setwd("/pathway/to/output/directory/")
# .rnk files from DoRothEA output by severity status
### Dendritic Cells files
Moderate_Rnk_FILE <- file.path(Folder,"Moderate_RnkFile.rnk")
Control_Rnk_FILE <- file.path(Folder,"Control_RnkFile.rnk")
Severe_Rnk_FILE <- file.path(Folder,"Severe_RnkFile.rnk")
### B-Cells files
B_Moderate_Rnk100_FILE <- file.path(Folder,"Moderate_Rnktop100.rnk")
B_Control_Rnk100_FILE <- file.path(Folder,"Control_Rnktop100.rnk")
B_Severe_Rnk100_FILE <- file.path(Folder,"Severe_Rnktop100.rnk")
# output directory for Webgestalt Results
outputDirectory <- getwd()
##################### DENDRITIC CELLS
############## MODERATE HEALTH STATUS
#Run Webgestalt with KEGG database in MODERATE health status
enrichResult_Moderate_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                      organism="hsapiens",
                                      enrichDatabase = "pathway_KEGG",
                                      enrichDatabaseType = "genesymbol",
                                      interestGeneFile= "Moderate_RnkFile.rnk",
                                      interestGeneType="genesymbol",
                                      isOutput=TRUE,
                                      sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                      gseaP = 1, setCoverNum = 10,
                                      outputDirectory=outputDirectory,
                                      projectName= "Moderate_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in MODERATE health status
enrichResult_Moderate_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                       organism="hsapiens",
                                       enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                       enrichDatabaseType = "genesymbol",
                                       interestGeneFile= "Moderate_RnkFile.rnk",
                                       interestGeneType="genesymbol",
                                       isOutput=TRUE,
                                       sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                       gseaP = 1, setCoverNum = 10,
                                       outputDirectory=outputDirectory,
                                       projectName= "Moderate_GSEA_GO-BP")

#Run Webgestalt with Reactome database in MODERATE health status 
enrichResult_Moderate_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                          organism="hsapiens",
                                          enrichDatabase = "pathway_Reactome",
                                          enrichDatabaseType = "genesymbol",
                                          interestGeneFile= "Moderate_RnkFile.rnk",
                                          interestGeneType="genesymbol",
                                          isOutput=TRUE,
                                          sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                          gseaP = 1, setCoverNum = 10,
                                          outputDirectory=outputDirectory,
                                          projectName= "Moderate_GSEA_Reactome")

#Run Webgestalt with Panther database in MODERATE health status  
enrichResult_Moderate_Panther <- WebGestaltR(enrichMethod="GSEA",
                                         organism="hsapiens",
                                         enrichDatabase = "pathway_Panther",
                                         enrichDatabaseType = "genesymbol",
                                         interestGeneFile= "Moderate_RnkFile.rnk",
                                         interestGeneType="genesymbol",
                                         isOutput=TRUE,
                                         sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                         gseaP = 1, setCoverNum = 10,
                                         outputDirectory=outputDirectory,
                                         projectName= "Moderate_RnkFile.rnk")

#Run Webgestalt with Wikipathway database in MODERATE health status
enrichResult_Moderate_Wikipathway <- WebGestaltR(enrichMethod="GSEA",
                                             organism="hsapiens",
                                             enrichDatabase = "pathway_Wikipathway",
                                             enrichDatabaseType = "genesymbol",
                                             interestGeneFile= "Moderate_RnkFile.rnk",
                                             interestGeneType="genesymbol",
                                             isOutput=TRUE,
                                             sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                             gseaP = 1, setCoverNum = 10,
                                             outputDirectory=outputDirectory,
                                             projectName= "Moderate_GSEA_Wikipathway")

############## CONTROL HEALTH STATUS
#Run Webgestalt with KEGG database in Control health status 
enrichResult_Crt_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                     organism="hsapiens",
                                     enrichDatabase = "pathway_KEGG",
                                     enrichDatabaseType = "genesymbol",
                                     interestGeneFile= "Control_RnkFile.rnk",
                                     interestGeneType="genesymbol",
                                     isOutput=TRUE,
                                     sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                     gseaP = 1, setCoverNum = 10,
                                     outputDirectory=outputDirectory,
                                     projectName= "Crt_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in Control health status 
enrichResult_Crt_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                      organism="hsapiens",
                                      enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                      enrichDatabaseType = "genesymbol",
                                      interestGeneFile= "Control_RnkFile.rnk",
                                      interestGeneType="genesymbol",
                                      isOutput=TRUE,
                                      sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                      gseaP = 1, setCoverNum = 10,
                                      outputDirectory=outputDirectory,
                                      projectName= "Crt_GSEA_GO-BP")

#Run Webgestalt with Reactome database in Control health status 
enrichResult_Crt_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                         organism="hsapiens",
                                         enrichDatabase = "pathway_Reactome",
                                         enrichDatabaseType = "genesymbol",
                                         interestGeneFile= "Control_RnkFile.rnk",
                                         interestGeneType="genesymbol",
                                         isOutput=TRUE,
                                         sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                         gseaP = 1, setCoverNum = 10,
                                         outputDirectory=outputDirectory,
                                         projectName= "Crt_GSEA_Reactome")

#Run Webgestalt with Panther database in Control health status
enrichResult_Crt_Panther <- WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase = "pathway_Panther",
                                        enrichDatabaseType = "genesymbol",
                                        interestGeneFile= "Control_RnkFile.rnk",
                                        interestGeneType="genesymbol",
                                        isOutput=TRUE,
                                        sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                        gseaP = 1, setCoverNum = 10,
                                        outputDirectory=outputDirectory,
                                        projectName= "Crt_GSEA_Panther")

#Run Webgestalt with Wikipathway database in Control health status
enrichResult_Crt_Wikipathway <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase = "pathway_Wikipathway",
                                            enrichDatabaseType = "genesymbol",
                                            interestGeneFile= "Control_RnkFile.rnk",
                                            interestGeneType="genesymbol",
                                            isOutput=TRUE,
                                            sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                            gseaP = 1, setCoverNum = 10,
                                            outputDirectory=outputDirectory,
                                            projectName= "Crt_GSEA_Wikipathway")

############## SEVERE HEALTH STATUS
#Run Webgestalt with KEGG database in Severe health status 
enrichResult_Svr_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                     organism="hsapiens",
                                     enrichDatabase = "pathway_KEGG",
                                     enrichDatabaseType = "genesymbol",
                                     interestGeneFile= "Severe_RnkFile.rnk",
                                     interestGeneType="genesymbol",
                                     isOutput=TRUE,
                                     sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                     gseaP = 1, setCoverNum = 10,
                                     outputDirectory=outputDirectory,
                                     projectName= "Svr_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in Severe health status
enrichResult_Svr_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                      organism="hsapiens",
                                      enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                      enrichDatabaseType = "genesymbol",
                                      interestGeneFile= "Severe_RnkFile.rnk",
                                      interestGeneType="genesymbol",
                                      isOutput=TRUE,
                                      sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                      gseaP = 1, setCoverNum = 10,
                                      outputDirectory=outputDirectory,
                                      projectName= "Svr_GSEA_GO-BP")

#Run Webgestalt with Reactome database in Severe health status 
enrichResult_Svr_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                         organism="hsapiens",
                                         enrichDatabase = "pathway_Reactome",
                                         enrichDatabaseType = "genesymbol",
                                         interestGeneFile= "Severe_RnkFile.rnk",
                                         interestGeneType="genesymbol",
                                         isOutput=TRUE,
                                         sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                         gseaP = 1, setCoverNum = 10,
                                         outputDirectory=outputDirectory,
                                         projectName= "Svr_GSEA_Reactome")

#Run Webgestalt with Panther database in Severe health status
enrichResult_Svr_Panther <- WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase = "pathway_Panther",
                                        enrichDatabaseType = "genesymbol",
                                        interestGeneFile= "Severe_RnkFile.rnk",
                                        interestGeneType="genesymbol",
                                        isOutput=TRUE,
                                        sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                        gseaP = 1, setCoverNum = 10,
                                        outputDirectory=outputDirectory,
                                        projectName= "Svr_GSEA_Panther")

#Run Webgestalt with Wikipathway database in Severe health status
enrichResult_Svr_Wikipathway <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase = "pathway_Wikipathway",
                                            enrichDatabaseType = "genesymbol",
                                            interestGeneFile= "Severe_RnkFile.rnk",
                                            interestGeneType="genesymbol",
                                            isOutput=TRUE,
                                            sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                            gseaP = 1, setCoverNum = 10,
                                            outputDirectory=outputDirectory,
                                            projectName= "Svr_GSEA_Wikipathway")


##################### B-CELLS
################### Moderate top 100
#Run Webgestalt with KEGG database in Moderate health status 
B_enrichResult_Moderate_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase = "pathway_KEGG",
                                            enrichDatabaseType = "genesymbol",
                                            interestGeneFile= "Moderate_Rnktop100.rnk",
                                            interestGeneType="genesymbol",
                                            isOutput=TRUE,
                                            sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                            gseaP = 1, setCoverNum = 10,
                                            outputDirectory=outputDirectory,
                                            projectName= "Moderate100_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in Moderate health status 
B_enrichResult_Moderate_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                             organism="hsapiens",
                                             enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                             enrichDatabaseType = "genesymbol",
                                             interestGeneFile= "Moderate_Rnktop100.rnk",
                                             interestGeneType="genesymbol",
                                             isOutput=TRUE,
                                             sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                             gseaP = 1, setCoverNum = 10,
                                             outputDirectory=outputDirectory,
                                             projectName= "Moderate100_GSEA_GO-BP")

#Run Webgestalt with Reactome database in Moderate health status  
B_enrichResult_Moderate_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                                organism="hsapiens",
                                                enrichDatabase = "pathway_Reactome",
                                                enrichDatabaseType = "genesymbol",
                                                interestGeneFile= "Moderate_Rnktop100.rnk",
                                                interestGeneType="genesymbol",
                                                isOutput=TRUE,
                                                sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                                gseaP = 1, setCoverNum = 10,
                                                outputDirectory=outputDirectory,
                                                projectName= "Moderate100_GSEA_Reactome")

#Run Webgestalt with Panther database in Moderate health status 
B_enrichResult_Moderate_Panther <- WebGestaltR(enrichMethod="GSEA",
                                               organism="hsapiens",
                                               enrichDatabase = "pathway_Panther",
                                               enrichDatabaseType = "genesymbol",
                                               interestGeneFile= "Moderate_Rnktop100.rnk",
                                               interestGeneType="genesymbol",
                                               isOutput=TRUE,
                                               sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                               gseaP = 1, setCoverNum = 10,
                                               outputDirectory=outputDirectory,
                                               projectName= "Moderate100_GSEA_Panther")

#Run Webgestalt with Wikipathway database in Moderate health status 
B_enrichResult_Moderate_Wikipathway <- WebGestaltR(enrichMethod="GSEA",
                                                   organism="hsapiens",
                                                   enrichDatabase = "pathway_Wikipathway",
                                                   enrichDatabaseType = "genesymbol",
                                                   interestGeneFile= "Moderate_Rnktop100.rnk",
                                                   interestGeneType="genesymbol",
                                                   isOutput=TRUE,
                                                   sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                                   gseaP = 1, setCoverNum = 10,
                                                   outputDirectory=outputDirectory,
                                                   projectName= "Moderate100_GSEA_Wikipathway")

################### Control top 100
#Run Webgestalt with KEGG database in Control health status  
B_enrichResult_Control_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                           organism="hsapiens",
                                           enrichDatabase = "pathway_KEGG",
                                           enrichDatabaseType = "genesymbol",
                                           interestGeneFile= "Control_Rnktop100.rnk",
                                           interestGeneType="genesymbol",
                                           isOutput=TRUE,
                                           sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                           gseaP = 1, setCoverNum = 10,
                                           outputDirectory=outputDirectory,
                                           projectName= "Control100_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in Control health status
B_enrichResult_Control_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                            enrichDatabaseType = "genesymbol",
                                            interestGeneFile= "Control_Rnktop100.rnk",
                                            interestGeneType="genesymbol",
                                            isOutput=TRUE,
                                            sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                            gseaP = 1, setCoverNum = 10,
                                            outputDirectory=outputDirectory,
                                            projectName= "Control100_GSEA_GO-BP")

#Run Webgestalt with Reactome database in Control health status 
B_enrichResult_Control_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                               organism="hsapiens",
                                               enrichDatabase = "pathway_Reactome",
                                               enrichDatabaseType = "genesymbol",
                                               interestGeneFile= "Control_Rnktop100.rnk",
                                               interestGeneType="genesymbol",
                                               isOutput=TRUE,
                                               sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                               gseaP = 1, setCoverNum = 10,
                                               outputDirectory=outputDirectory,
                                               projectName= "Control100_GSEA_Reactome")

#Run Webgestalt with Panther database in Control health status
B_enrichResult_Control_Panther <- WebGestaltR(enrichMethod="GSEA",
                                              organism="hsapiens",
                                              enrichDatabase = "pathway_Panther",
                                              enrichDatabaseType = "genesymbol",
                                              interestGeneFile= "Control_Rnktop100.rnk",
                                              interestGeneType="genesymbol",
                                              isOutput=TRUE,
                                              sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                              gseaP = 1, setCoverNum = 10,
                                              outputDirectory=outputDirectory,
                                              projectName= "Control100_GSEA_Panther")

#Run Webgestalt with Wikipathway database in Control health status 
B_enrichResult_Control_Wikpathway <- WebGestaltR(enrichMethod="GSEA",
                                                 organism="hsapiens",
                                                 enrichDatabase = "pathway_Wikipathway",
                                                 enrichDatabaseType = "genesymbol",
                                                 interestGeneFile= "Control_Rnktop100.rnk",
                                                 interestGeneType="genesymbol",
                                                 isOutput=TRUE,
                                                 sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                                 gseaP = 1, setCoverNum = 10,
                                                 outputDirectory=outputDirectory,
                                                 projectName= "Control100_GSEA_Wikipathway")

################### Severe top 100
#Run Webgestalt with KEGG database in Severe health status 
B_enrichResult_Severe_KEGG <- WebGestaltR(enrichMethod="GSEA",
                                          organism="hsapiens",
                                          enrichDatabase = "pathway_KEGG",
                                          enrichDatabaseType = "genesymbol",
                                          interestGeneFile= "Severe_Rnktop100.rnk",
                                          interestGeneType="genesymbol",
                                          isOutput=TRUE,
                                          sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                          gseaP = 1, setCoverNum = 10,
                                          outputDirectory=outputDirectory,
                                          projectName= "Severe100_GSEA_KEGG")

#Run Webgestalt with GO-BP noRedundant database in Severe health status 
B_enrichResult_Severe_GO-BP <- WebGestaltR(enrichMethod="GSEA",
                                           organism="hsapiens",
                                           enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                           enrichDatabaseType = "genesymbol",
                                           interestGeneFile= "Severe_Rnktop100.rnk",
                                           interestGeneType="genesymbol",
                                           isOutput=TRUE,
                                           sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                           gseaP = 1, setCoverNum = 10,
                                           outputDirectory=outputDirectory,
                                           projectName= "Severe100_GSEA_GO-BP")

#Run Webgestalt with Reactome database in Severe health status
B_enrichResult_Severe_Reactome <- WebGestaltR(enrichMethod="GSEA",
                                              organism="hsapiens",
                                              enrichDatabase = "pathway_Reactome",
                                              enrichDatabaseType = "genesymbol",
                                              interestGeneFile= "Severe_Rnktop100.rnk",
                                              interestGeneType="genesymbol",
                                              isOutput=TRUE,
                                              sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                              gseaP = 1, setCoverNum = 10,
                                              outputDirectory=outputDirectory,
                                              projectName= "Severe100_GSEA_Reactome")

#Run Webgestalt with Panther database in Severe health status 
B_enrichResult_Severe_Panther <- WebGestaltR(enrichMethod="GSEA",
                                             organism="hsapiens",
                                             enrichDatabase = "pathway_Panther",
                                             enrichDatabaseType = "genesymbol",
                                             interestGeneFile= "Severe_Rnktop100.rnk",
                                             interestGeneType="genesymbol",
                                             isOutput=TRUE,
                                             sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                             gseaP = 1, setCoverNum = 10,
                                             outputDirectory=outputDirectory,
                                             projectName= "Severe100_GSEA_Panther")

#Run Webgestalt with Wikipathway database in Severe health status 
B_enrichResult_Severe_Wikipathway <- WebGestaltR(enrichMethod="GSEA",
                                                 organism="hsapiens",
                                                 enrichDatabase = "pathway_Wikipathway",
                                                 enrichDatabaseType = "genesymbol",
                                                 interestGeneFile= "Severe_Rnktop100.rnk",
                                                 interestGeneType="genesymbol",
                                                 isOutput=TRUE,
                                                 sigMethod="top", topThr=10, minNum=5, maxNum = 500, reportNum = 30, perNum = 10000,
                                                 gseaP = 1, setCoverNum = 10,
                                                 outputDirectory=outputDirectory,
                                                 projectName= "Severe100_GSEA_Wikipathway")

############ PLOT THE OVERALL Results
# By manually exploring the enrichment tables we filtered the results by selecting the pathways with a biological
# relationship in COVID-19, Here we add the resulting tables. However, the user can also manually explore the whole results
# obtained above.
# Bar plot de los datos de GSEA
############ DENDRITIC CELLS
# Library for the enrichmentplots
library(ggplot2)
# MODERATE GSEA Results
GSEA_Mld = read.delim("GSEA_mld.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
#REMOVE pathways with no apparent biological relevance
GSEA_Mld <- GSEA_Mld[-c(8, 13, 14), ]
#CONTROL GSEA Results
GSEA_Crt = read.delim("GSEA_Ctr.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
#REMOVE pathways with no apparent biological relevance
GSEA_Crt <- GSEA_Crt[-c(2, 7, 11, 12,14,15,16,20,22,23,24,27,28,29,30,32,34,36), ]
#SEVERE GSEA Results
GSEA_Svr = read.delim("GSEA_Svr.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
#REMOVE pathways with no apparent biological relevance
GSEA_Svr <- GSEA_Svr[-c(1, 11, 12, 13, 23,24,25,27,28,29,30), ]

#Enrichment plot by activation or inactivation
GSEA_Mld$Status <- ifelse(GSEA_Mld$NES >=0, "Active", "Inactive")
GSEA_Crt$Status <- ifelse(GSEA_Crt$NES >=0, "Active", "Inactive")
GSEA_Svr$Status <- ifelse(GSEA_Svr$NES >=0, "Active", "Inactive")
#MODERATE
S1 <- ggplot(GSEA_Mld, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database, x=-1), color="white", fontface="bold") +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S1
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S1<- S1 + theme(axis.title.y = element_blank())
S1
#CONTROL
S2 <- ggplot(GSEA_Crt, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database, x=-1), color="white", fontface="bold") +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S2
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S2<- S2 + theme(axis.title.y = element_blank())
S2
#SEVERE
S3 <- ggplot(GSEA_Svr, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database, x=1), color="white", fontface="bold") +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S3
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S3<- S3 + theme(axis.title.y = element_blank())
S3
## REMOVE LEGENDS
S2_NL <- S2 + theme(legend.position = "none")
S2_NL
#FINAL ENRICHMENT PLOT
left_col_COLs<-cowplot::plot_grid(S1,labels = c("A"),ncol = 1,rel_heights = c(0.5,0.5))
top_panel_NL_COLs<-cowplot::plot_grid(left_col_COLs,S2_NL,labels = c("","B"),ncol = 1,rel_heights = c(1,1.5))
finalplot_COLs<-cowplot::plot_grid(top_panel_NL_COLs,S3,labels = c("","C"),ncol = 1,rel_heights = c(1,0.5))
finalplot_COLs
ggsave(finalplot_COLs,file="plotGSEA_COL.png",width = 20,height =25,units = "cm",dpi = 300)

############## Enrichment plot for B-cells
# CONTROL GSEA Results
GSEA_C_Bcells = read.delim("GSEA_Control_BCells.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
# SEVERE GSEA Results
GSEA_S_Bcells = read.delim("GSEA_Severe_BCells.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
# MODERATE GSEA Results
GSEA_M_Bcells = read.delim("GSEA_Moderate_BCells.txt", header = T, stringsAsFactors = FALSE, quote = "", sep = "\t")
#Enrichment plot by activation or inactivation
GSEA_C_Bcells$Status <- ifelse(GSEA_C_Bcells$NES >=0, "Active", "Inactive")
GSEA_C_Bcells
GSEA_S_Bcells$Status <- ifelse(GSEA_S_Bcells$NES >=0, "Active", "Inactive")
GSEA_S_Bcells
GSEA_M_Bcells$Status <- ifelse(GSEA_M_Bcells$NES >=0, "Active", "Inactive")
GSEA_M_Bcells
#MODERATE
S4 <- ggplot(GSEA_M_Bcells, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database), color="white", fontface="bold",position = position_stack(vjust = 0.5)) +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S4
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S4<- S4 + theme(axis.title.y = element_blank())
S4
#CONTROL
S5 <- ggplot(GSEA_C_Bcells, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database, x=1), color="white", fontface="bold") +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S5
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S5<- S5 + theme(axis.title.y = element_blank())
S5
# NO legend
S5_NL <- S5 + theme(legend.position = "none")
S5_NL
#SEVERE
S6 <- ggplot(GSEA_S_Bcells, aes(x= NES, y= reorder(Pathway, p.adjust), fill=Status)) +
  geom_bar(stat = "identity")  +
  xlab("NES") + ylab("Pathway") +
  scale_fill_manual("Status", values = c("Active" = "#FF2600", "Inactive" = "#0433FF")) +
  geom_text(aes(label=database), color="white", fontface="bold", position = position_stack(vjust = 0.5)) +
  theme(legend.direction = "vertical", legend.box = "vertical") +
  theme_bw()
S6
# Key function: use element_blank() to suppress axis labels. (x label: axis.title.x = element_blank())
S6<- S6 + theme(axis.title.y = element_blank())
S6
# NO legend
S6_NL <- S6 + theme(legend.position = "none")
S6_NL
#FINAL ENRICHMENT PLOT
left_col_COLs_Bcells<-cowplot::plot_grid(S4,labels = c("A"),ncol = 1,rel_heights = c(1,1))
top_panel_NL_COLs_Bcells<-cowplot::plot_grid(left_col_COLs_Bcells,S5_NL,labels = c("","B"),ncol = 1,rel_heights = c(1,0.8))
finalplot_COLs_Bcells<-cowplot::plot_grid(top_panel_NL_COLs_Bcells,S6_NL,labels = c("","C"),ncol = 1,rel_heights = c(1,1))
finalplot_COLs_Bcells
ggsave(finalplot_COLs_Bcells,file="plotGSEA_COL_Bcells.png",width = 20,height =25,units = "cm",dpi = 300)


