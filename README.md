# Covid-19_scRNAseq
 In this section we present an analysis based on data of single cell RNAseq coming from patients with COVID-19 stratified through 3 levels or clinical responses: normal, covid-19 moderate, and covid-19 severe or critical. Data was obtained from this paper https://www.nature.com/articles/s41591-020-0901-9. In this paper the authors proceed with a BronchoALveolar lavage Fluid (BALF) of 12 patients (3 moderate M1-M3, 6 patients wit severe/critical infection S1-S6, and 3 normal controls-HC1-HC3).

# XGBoost tree classification.

 Our machine learning analysis start from the count matrix, integrating all the count of genes identified in all the 90,000 unique cells. The main files are:
 1) Data_filtered.txt,Data_filtered.txt ( expression matrix as "sparse Matrix")
 2) Labels.txt (label of the samples ordered as indicated in the count matrix).
 3) genes.txt ( list of genes orderedas the count matrix).
 
Requirements: Install 
```
from scipy import sparse
from numpy import array
import numpy as np
import pandas as pd
from numpy import savetxt 
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.externals import joblib
print('Libraries Imported')
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.externals import joblib
print('Libraries Imported')
from scipy import sparse
from numpy import array
from numpy import savetxt 
import pickle
import xgboost as xgb
import matplotlib.pylab as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix
from matplotlib import pyplot
import matplotlib
```

Let's star the analysis.

# Files preprocessing.

Given that matrix is written in sparse notation, we separate the varibles to build the normal matrix. In sparse notation the first number into a row is the row of the normal matrix, the second number on the row is the colums; and the third is the numerical values of the entry. We separate and define files to these information. It was done in bash to define these files:
1) rows.csv (File with the rows)
2) cols.csv (File with columns)
3) values.csv (File with the values at row and columns specified in rows.csv and cols.csv files)

```
osbaldo@nautilus:~/COVID-19/singleCell$ sed -n '3,143014366p' Data_filtered.txt | awk '{print$1}' > rows.csv
osbaldo@nautilus:~/COVID-19/singleCell$ sed -n '3,143014366p' Data_filtered.txt | awk '{print$2}' > cols.csv
osbaldo@nautilus:~/COVID-19/singleCell$ sed -n '3,143014366p' Data_filtered.txt | awk '{print$3}' > values.csv
```

We star the analysis with these files:
1) cols.csv  
2) Data_filtered.txt  
3) genes.tsv  
4) labels.tsv  
5) medical_class.csv  (This file contains the classification of the patients: 0 Normal, 1 COVID-modest; 2) COVID-19 severe or response)
6) rows.csv  
7) values.csv


Run this script in python to buld the model for machine learinig analysis 
```
python3 model.py
```
The genrated files are located in the folders: 
1) /dat/  (contain the X and y picled)
2) /Splited_data/ (contain X_test.pickle  X_train.pickle  y_test.pickle  y_train.pickle)

