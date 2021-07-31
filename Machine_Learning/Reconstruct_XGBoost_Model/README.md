# Generation of XGBoost Model  
This section, we present the procedure by which we construct the machine learning model to classify cells of patients with COVID-19 with a moderate or severe response from single-cell RNASeq samples. From the original paper the samples were stratified through 2 levels or clinical responses: covid-19 moderate (0), and covid-19 severe or critical (1). We start with these files

## Initial files
|Files  |      Description  | 
|------------|:---------------:|
|  matrix.zip          |  Count matrix of the cells         | 
|  genes.tsv          |  Lists of genes        | 
|  barcodes.zip          | Classification of the cells as moderate (0) or severe(1) responses in COVID-19 patients       |
| col.zip | Index for columns  |
| rows.zip        |  Index for rows                   |
| values.zip| Index for numric values of count matrix |
 

## Scripts for generating the model
In addition we have four scripts in python, these required to build and test the model

|Files  |      Description  | 
|------------|:---------------:|
|  model.py          |  Required to format the data for training the model         | 
|  load_libraries.py          |  Concomitant file to load the libraries in python        |
|mlcovid.py| Traing the XGBoost model|
|crossvalidation.py| Crossvalidation test|

To build the XGBoost model we proceeded as follows:

**Step 1** Run in command line the script:
```
> python3 model.py
```
We will obtain the data in a proper notation and these separated in the independent (X) and depende variables (y). These variables are compressed in the file X_y.zip. A similar code was a applied to format the validation data, these obtained from [Ren et al](https://www.cell.com/cell/fulltext/S0092-8674(21)00148-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421001483%3Fshowall%3Dtrue). These variables are included int X_val_y_val.zip
 
**Step 2** Open a sesion of python and load all the libraries required. In terminal and run:

```
> python3 mlcovid.py
```

### Step 3

Open a sesion of python and load all the libraries required. In terminal and run:
```
> python3 crossvalidation.py
```

## Outputs


## Figures



