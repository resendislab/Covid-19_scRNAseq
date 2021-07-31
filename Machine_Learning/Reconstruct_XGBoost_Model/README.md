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

### Step 1

Run in command line the script:

```
> python3 model.py
```
We will obtain the data in a proper notation and these separated in the independent (X) and depende variables (y).
 
### Step 2.

Open a sesion of python and load all the libraries required. In terminal and run:

```
> python3 mlcovid.py
```

### Step 3

Open a sesion of python and load all the libraries required. In terminal and run:
```
> python3 crossvalidation.py
```





