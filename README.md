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
