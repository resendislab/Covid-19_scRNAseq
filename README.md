# Single Cell RNA-Seq Data Analysis for COVID-19 samples

## Attribution
The pipeline was developed in the [Human Systems Biology Group](https://resendislab.github.io/) at [INMEGEN](https://www.inmegen.gob.mx/) that you can [cite].

## Getting Started
The pipeline was used to 1. develop a severity classification method for COVID-19 patients using machine learning techniques. 2. study dysregulated pathways in COVID-19 lung cells. If you have any questions contact us at avazquezj(at)inmegen.gob.mx and oresendis(at)inmegen.gob.mx.

The full pipeline is structured in the three steps:

* Data Download\
Raw scRNA-seq data of Covid-19 patients were obtained from GEO (GSE145926). To sum up, bronchoalveolar lavage fluid (BALF) cells were collected from 12 patients and grouped according to their symptoms as healthy, moderate, and severe. Samples were sequenced using 10x Genomics technology [(Liao et al., 2020)](https://www.nature.com/articles/s41591-020-0901-9).\
The associated scrips can be found in the Download_Data folder:
    - The Download_data.R file downloads data and save them in a local folder.    
    - The Load_data.R file merges the data with labels for control, moderate, and severe files into one sparse matrix.
    
    For the machine learning analysis, three files are saved: the matrix count, gene names, and label samples. In the case of the immune landscape analysis, a Seurat object is saved containing the data and their metadata.

* [Machine Learning](Machine_Learning/README.md)\
To improve patients' prognosis and treatment, we proposed a classification tree using 11 genes taken from moderate and severe patients. This analysis was done in python.


* [Immune landscape analysis](Immune_Landscape/README.md)\
We explored the differences in the immune response among patients. We used R to conduct singlge cell analysis.
