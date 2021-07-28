
# XGBoost Model

In the paper we have constructed a XGBoost model with capacities to classify cell from moderate and severe COVID-19 patients from their scRNASeq profiles. All the high-throughput data were obtained from the Liao et al [paper](https://www.nature.com/articles/s41591-020-0901-9). With the purpose to assess the capacities of this machine learning techniques, we have constructed two models. One of the models have taked into account the entire set of genes that were reported in the original paper. The second model has reduced the set genes by excluding those genes associated with quality factors. The material is organized by two folders, one specifying the completed and reduced model respectively. Combining both folders with the scrip Fig2.py we obtained this figure
![confuse](Fig2.png) 
