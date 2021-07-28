
# Load validation data. Check that you are in a folder with load_libraries
from load_libraries import *
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

def myfun(self=None):
    return model_bytearray

pickle_in = open("model_reduced.pickle","rb")
model_reduced = pickle.load(pickle_in);

#Load the files test
pickle_in = open("X_test.pickle","rb")
X_test = pickle.load(pickle_in);

pickle_in = open("y_test.pickle","rb")
y_test = pickle.load(pickle_in);


pickle_in = open("X_val.pickle","rb")
X_val = pickle.load(pickle_in);

pickle_in = open("y_val.pickle","rb")
y_val = pickle.load(pickle_in);


# Write in D format
D_test = xgb.DMatrix(X_test, label=y_test);
D_val = xgb.DMatrix(X_val, label=y_val);


# Performance of the model evaluating the test dataset.
pred_test= model_reduced.predict(D_test)
pred_val= model_reduced.predict(D_val)

print(classification_report(y_test, pred_test))
cm_test = confusion_matrix(y_test, pred_test)
cm_val = confusion_matrix(y_val, pred_val)

## To generate ONLY the confusion matrix with labels
cmap = plt.get_cmap('Blues')

cm_display = ConfusionMatrixDisplay(confusion_matrix=cm_test,display_labels=['Moderate','Severe']).plot(cmap = cmap,colorbar=False)
plt.savefig('confusematrix_test.png',dpi=800)

cm_display = ConfusionMatrixDisplay(confusion_matrix=cm_val,display_labels=['Moderate','Severe']).plot(cmap = cmap,colorbar=False)
plt.savefig('confusematrix_val.png',dpi=800)

