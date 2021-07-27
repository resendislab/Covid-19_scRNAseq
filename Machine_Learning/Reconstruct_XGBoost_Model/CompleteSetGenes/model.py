# STEP: LOAD LIBRARIES
print('Processing...')
from scipy import sparse
from numpy import array
import numpy as np
import pandas as pd
from numpy import savetxt 
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
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
import os
print('Libraries Imported')
#######################################################################################
# Preparation of Matrix 
r = pd.read_csv('rows.csv', header = None);
c = pd.read_csv('cols.csv',header = None);
v = pd.read_csv('values.csv',header = None);
colname = pd.read_csv('genes.tsv', header= None);
rowname = pd.read_csv('barcodes.tsv',header= None);
r1 = r.values.flatten()
c1= c.values.flatten()
v1 = v.values.flatten()
colname = colname.values.flatten()
rowname = rowname.values.flatten()
#Dense Matrix
b= sparse.coo_matrix((v1,(r1,c1)));
dense = b.toarray();
dense.shape
z = pd.DataFrame(data = dense[1:,1:],index=rowname,columns=colname);
data = z

#####################################################################################
# Data preparation
y = data['Severity'];
X= data.drop(['Severity'],axis=1);
X.index = range(len(X))
y.index = range(len(y))

#####################################################################################
## Save files as pickle
os.system('mkdir dat')
pickle.dump(X, open("dat/X.pickle", 'wb'), protocol=4);
pickle_out = open("dat/y.pickle","wb")
pickle.dump(y, pickle_out)
pickle_out.close()

######################################################################################
print('Finish')

