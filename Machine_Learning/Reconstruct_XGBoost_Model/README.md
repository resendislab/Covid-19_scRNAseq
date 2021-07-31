# Generation of XGBoost Model  
This section, we present the procedure by which we construct the machine learning model to classify cells of patients with COVID-19 with a moderate or severe response from single-cell RNASeq samples. From the original paper the samples were stratified through 2 levels or clinical responses: covid-19 moderate (0), and covid-19 severe or critical (1). We start with these files

|Files  |      Description  | 
|------------|:---------------:|
|  matrix.zip          |  Count matrix of the cells         | 
|  genes.tsv          |  Lists of genes        | 
|  barcodes.zip          | Classification of the cells as moderate (0) or severe(1) responses in COVID-19 patients       |
| col.zip | Index for columns  |
| rows.zip        |  Index for rows                   |
| values.zip| Index for numric values of count matrix |
 


In addition we have four scripts in python, these required to build and test the model

|Files  |      Description  | 
|------------|:---------------:|
|  model.py          |  Required to format the data for training the model         | 
|  load_libraries.py          |  Concomitant file to load the libraries in python        |
|mlcovid.py| Traing the XGBoost model|
|crossvalidation.py| Crossvalidation test|

To build the XGBoost model we proceeded as follows:

## Step 1

Run in command line the script:

```
> python3 model.py
```
We will obtain the data in a proper notation and these separated in the independent (X) and depende variables (y).
 
## Step 2.

Open a sesion of python and load all the libraries required. In terminal and run:

```
> python3 mlcovid.py
```

## Step 3

Open a sesion of python and load all the libraries required. In terminal and run:
```
> python3 crossvalidation.py
```


## Explaning analysis

In this setion we explore the set of genes that played a significant role in differentiate the samples at any level (normal, moderate and severe COVID-19). The first plot render the genes that have more contribution to separate the samples. 

![confuse](confusematrixblue.png)
```
import os
import pickle
import shap

```
Lets start th SHAP analysis
```
import os
# SHAP analysis
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)
os.system('mkdir pickle')
# Save shap_values in pickle
pickle.dump(shap_values, open("pickle/shap_values.pickle", 'wb'), protocol=4);

#load shap_values. This will allow to load the shao_values file withour run again the code.
import os
os.system('mkdir pickle')
pickle_in = open("pickle/shap_values.pickle","rb")
shap_values = pickle.load(pickle_in);
```
Finally we recap and save the figure

```

shap.summary_plot(shap_values, X)
pyplot.savefig('figures/shap_summary1.png')
```
In figure below we identify those genes with highest relevance in classify the patients as normal (0), moderate (1), and severe (2) COVID-19 phenotype. As Figure shows, genes have different degree of contribution in most of the subgroups. Besides, only a few genes have a specific contitribution in some groups. 

![summary](shap_summary1.png)


To have a more detail description of the role of important genes into the clasification, we procedded to analyze the profile behavior at each classitication for each important gene shown in the previous plot.  

```
##### SUMMANRY PLOTS
import matplotlib.pyplot as plt
for i in range(3):
    pl.clf()
    shap.summary_plot(shap_values[i], X)
    plt.savefig('sum_' + str(i)+ '.png', format='png', dpi=300, bbox_inches='tight')

```
For instance, in the case of the classification for patients with COVID-19 severe, the gene expression profile of the genes are given below:

![dependence](/summary/sum_2.png)

## Dependence plot

The last part of this analysis is focused to identify association between variables. We depicted these plots, called dependence plots, for main variables involved in classification. Below, we show one of this figures, but the entire set are stored in folder Dependence_figures.  
To obtain the dependence plots: 

```
## Figures
## We have selected those genes with importance in the classification.
lista = [12962,20721,24024,16163,16164,4903,25601,25602,8242,24024,12964,19123,25603,11543,25605,24762,13668,24587,2576,6044,20797]

os.chdir('figures')
for name in X_train.columns[lista]:
    shap.dependence_plot(name, shap_values[0], X, display_features=X_display)
    plt.savefig(str(name) + '.png', format='png', dpi=300, bbox_inches='tight')
```

![dependence](/Dependence_figures/TAOK1.png)



