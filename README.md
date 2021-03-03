# Single Cell RNA-Seq Data Analysis for COVID-19 samples

## Attribution
The pipeline was developed in the [Human Systems Biology Group](https://resendislab.github.io/) at [INMEGEN](https://www.inmegen.gob.mx/) that you can [cite].

## Getting Started
The pipeline was used to: 1. develop a severity classifcation method for COVID-19 patients usign machine learning techniques. 2. study disregulated pathways in COVID-19 lungs cells. If you have any questions contact us at avazquezj(at)inmegen.gob.mx and oresendis(at)inmegen.gob.mx.

The full pipeline is structured in the three steps:

* Data Download\
Raw scRNA-seq data of Covid-19 patients were obtained from GEO (GSE145926). To sum up, bronchoalveolar lavage fluid (BALF) cells were collected from 12 patients and grouped according to their symptoms as healthy, moderate, and severe. Samples were sequenced using 10x Genomics technology [(Liao et al., 2020)](https://www.nature.com/articles/s41591-020-0901-9).\
The associated scrips can be found in Download_Data folder:
    - The Download_data.R file downloads data and save them in a local folder.    
    - The Load_data.R file merge the data with labels for normal, moderate and severe files are associated into one sparse matrix.
    
    For the machine learning analysis three files are saved: the matrix count, gene names and labels samples. In the case for the immune landscape analysis, a seurat object is saved containing the data and the samples metadata.

* [Machine Learning](Machine_Learning/README.md)\
To improve patients prognosis and improve treatment, we proposed a classification tree using 11 genes taken from moderate and severe patients.


* [Immune landscape analysis](Immune_Landscape/README.md)\
We explored the differences in the immune response among patients.
