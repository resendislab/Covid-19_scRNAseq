#############################Crossvalidation################################################################################
############################################################################################################################
############################################################################################################################

# Import libraries

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import confusion_matrix, roc_auc_score
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import seaborn as sns
import shap
from load_libraries import *
import os

###########################################################################################################################
# Load files X and y.
print("Loading X and y")
pickle_in = open("dat/X.pickle","rb")
X = pickle.load(pickle_in);
pickle_in = open("dat/y.pickle","rb")
y = pickle.load(pickle_in);

###########################################################################################################################
# Model reconstruction
# seed to reproduce the data
np.random.seed(49)
# n_splits = 5
kfolds = KFold(n_splits=5, shuffle=True)

# Parameters
params = {
    'max_depth': 4,
    'eta': 0.2,
    'eval_metric':'auc',
    'objective': 'binary:hinge',
    'n_gpus': 0
}

steps = 1000

###########################################################################################################################
###STEP Crossvalidation
cnf =[] # list()
auc = [] # list()
features = X.columns.values.tolist()
thres = 0.5

count=0
os.system('mkdir crossvalidation')
for train_idx, test_idx in kfolds.split(X):
    count = count + 1;
    print(count)
    X_train, y_train = X.iloc[train_idx], y.iloc[train_idx]
    X_test, y_test = X.iloc[test_idx], y.iloc[test_idx]
    xg_train = xgb.DMatrix(X_train.values, feature_names= features, label=y_train.values)
    xg_test = xgb.DMatrix(X_test.values, feature_names=features, label=y_test.values)
    watchlist = [(xg_train, 'train'), (xg_test, 'test')]
    params = {
        'max_depth': 4,
        'eta': 0.2,
        'eval_metric':'auc',
        'objective': 'binary:hinge', 
        'seed':count,
        'n_gpus': 0
    }
    steps = 1000;
    print("Constructing model...")
    bst = xgb.train(params, xg_train, steps, watchlist, verbose_eval=False)
    preds = bst.predict(xg_test)
    a = confusion_matrix(y_test, (preds > thres).astype(int))
    b = roc_auc_score(y_test, preds)
    auc.append(b)
    del a
    del b 


# Array the model.
    bst_bytearray = bst.save_raw()[4:]
    def myfun(self=None):
        return bst_bytearray

    bst.save_raw = myfun


#####################################################################################################################################
# STEP. This is to analyze all the crossvalidation plots and see the variation in the list of genes
    explainer = shap.TreeExplainer(bst)
    shap_values = explainer.shap_values(X)
    vals = np.abs(shap_values).mean(0); # values importance list
    fi = pd.DataFrame(list(zip(X_train.columns,vals)),columns=['colname','feature_importance']);
    fi.sort_values(by=['feature_importance'],ascending=False,inplace=True);
    f1=fi.head(20);
    f1.to_csv('crossvalidation/shap_summary_' + str(count) + '.csv', sep =" ",header= True,index=False);
    X_display = X 
    plt.clf()
    shap.summary_plot(shap_values, X)
    pyplot.savefig('.../crossvalidation/shap_summary_' + str(count) + '.png',format='png', dpi=300, bbox_inches='tight')




######################################################################################################################################
# Save the data
AUC =pd.DataFrame(auc)
AUC.to_csv('crossvalidation/auc.csv')

plt.clf()
AUC.plot()
pyplot.savefig('crossvalidation/plot_AUC'  + '.png',format='png', dpi=300, bbox_inches='tight')


plt.clf()
AUC.hist(bins=20)
pyplot.savefig('crossvalidation/hist_AUC'  + '.png',format='png', dpi=300, bbox_inches='tight')

