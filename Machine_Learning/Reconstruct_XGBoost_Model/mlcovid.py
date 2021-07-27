
################################################################################################################
# Load X(data) and y(classes) files.
print('Libraries')
from load_libraries import *

# Load data (X and y)
pickle_in = open("...dat/X.pickle","rb")
X = pickle.load(pickle_in);
pickle_in = open("...dat/y.pickle","rb")
y = pickle.load(pickle_in);

# Load validation data.
pickle_in = open(".../X_val.pickle","rb")
X_val = pickle.load(pickle_in);
pickle_in = open(".../y_val.pickle","rb")
y_val = pickle.load(pickle_in);

###############################################################################################################
#Split the data and Save X, y train and test.
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state = 43,stratify = y);

import pickle
import os

os.system('mkdir Splited_data')
pickle.dump(X_train, open("Splited_data/X_train.pickle", 'wb'), protocol=4);
pickle.dump(X_test, open("Splited_data/X_test.pickle", 'wb'), protocol=4);

pickle_out = open("Splited_data/y_train.pickle","wb")
pickle.dump(y_train, pickle_out)
pickle_out.close()

pickle_out = open("Splited_data/y_test.pickle","wb")
pickle.dump(y_test, pickle_out)
pickle_out.close()

################################################################################################################
# Construct the XGBoost model. 
print('Constructing model..')

D_train = xgb.DMatrix(X_train, label=y_train);
D_test = xgb.DMatrix(X_test, label=y_test);

# Notation to validation test.
D_val = xgb.DMatrix(X_val, label=y_val);

#parameters
params = {
    'max_depth': 4,
    'eta': 0.2,
    'eval_metric':'auc',
    'objective': 'binary:hinge',  
    'seed':90,
    'n_gpus': 0
}

steps = 1000

model = xgb.train(params, D_train, steps)
#Prepare the format
model_bytearray = model.save_raw()[4:]
def myfun(self=None):
    return model_bytearray

model.save_raw = myfun

##################################################################################################################
# Save model
os.system('mkdir Model')
pickle_out = open("Model/model.pickle","wb")
pickle.dump(model, pickle_out)
pickle_out.close()

###################################################################################################################
# Assessment of the model. Confusion matrix
import os
import shap

# Performance of the model evaluating the test dataset.
pred= model.predict(D_test)
print(classification_report(y_test, pred))
cm = confusion_matrix(y_test, pred)
cm

# Performance of the model evaluating the VALIDATION dataset.
pred_val= model.predict(D_val)
print(classification_report(y_val, pred_val))
cm_val = confusion_matrix(y_val, pred_val)
cm_val



cmap = plt.get_cmap('Blues')

def plot_confusion_matrix(cm, classes, normalized=True, cmap=cmap): # 'bone'
    plt.figure(figsize=[7, 6])
    norm_cm = cm
    if normalized:
        norm_cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        sns.heatmap(norm_cm, annot=cm, fmt='g', xticklabels=classes, yticklabels=classes, cmap=cmap)

import os
os.system('mkdir figures')

# I will save the figure data to make the figure in the paper.
# training and test performance
pickle_out = open("figures/cm.pickle","wb")
pickle.dump(cm, pickle_out)
pickle_out.close()

plot_confusion_matrix(cm, ['Cov_Moderate', 'Cov-severe'])
pyplot.savefig('figures/confusematrix.png',dpi=800) 
print('Confuse matrix saved!')

# validation performance
pickle_out = open("figures/cm_val.pickle","wb")
pickle.dump(cm_val, pickle_out)
pickle_out.close()

plot_confusion_matrix(cm_val, ['Cov_Moderate', 'Cov-severe'])
pyplot.savefig('figures/confusematrix_validation.png',dpi=800) 
print('Confuse matrix saved!')


####################################################################################################################
# Figures of the paper.

import matplotlib.pylab as plt
import os
import pickle
import shap

# SHAP analysis
print('Shap plots in progress..')
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X) 


# Save shap_values
os.system('mkdir pickle_shap')
pickle.dump(shap_values, open("pickle_shap/shap_values.pickle", 'wb'), protocol=4);

X_display = X

plt.clf()
shap.summary_plot(shap_values, X)
pyplot.savefig('figures/shap_summary1.png',format='png', dpi=800, bbox_inches='tight')

plt.clf()
shap.summary_plot(shap_values, X, plot_type="bar")
pyplot.savefig('figures/shap_summary1_bar.png',format='png', dpi=800, bbox_inches='tight')
print('Shap plots saved!!')

######################################################################################################################################### 


