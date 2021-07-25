
# XGBoost Model

 In this section, we present the machine learning analysis of single-cell RNASeq samples for patients with COVID-19. From the original paper the samples were stratified through 2 levels or clinical responses: covid-19 moderate, and covid-19 severe or critical.
 
# Analysis severe(1) vs moderate(0) response. 

In this section we compared and identify genes that separate the behavior of Severe and moderate patients with COVID-19. We started the analysis from seven files: (server: /media/usb/osbaldo/COVID-19/singlecell/comparison_moderate_severe_covid_19)

    1) Data_filtered.txt (original de single cell).
    2) genes.tsv         (original de single cell).
    3) labels.tsv        (original de single cell).
    From these files we obtain:
    1) cols.csv
    5) medical_class.csv (This file contains the classification of the patients: 0 Normal, 1 COVID-moderate; 2) COVID-19 severe response)
    6) rows.csv
    7) values.csv

    In addition we have these scripts in python:
    1) model.py
    2) load_libraries.py
    3) mlcovid.py
    4) crossvalidation.py


In this situation the number o samples to test are:

|Tables Test |      Phenotype  | Samples |
|------------|:---------------:|--------:|
|            |  severe         |  12699  |
|            |  moderate       |   7014  |


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

## Step 1

Open a sesion of python and load all the libraries required. In terminal and run:
```
> python3 crossvalidation.py
```

