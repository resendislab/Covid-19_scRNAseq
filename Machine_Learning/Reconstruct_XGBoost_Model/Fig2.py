
######################################################################################################
######################################################################################################
######################################################################################################
# import libraries
import pickle
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import shap
from venn import venn
import matplotlib.pylab as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

# Load the files
pickle_in = open(".../X.pickle","rb")
X = pickle.load(pickle_in);

pickle_in = open(".../pickle_shap/shap_values.pickle","rb")
shap_values = pickle.load(pickle_in);

pickle_in = open(".../figures/cm.pickle","rb")
cm = pickle.load(pickle_in);

pickle_in = open(".../figures/cm_val.pickle","rb")
cm_val = pickle.load(pickle_in);


cmap = plt.get_cmap('Blues')



lista ={"Realization_1":{"RPLP1","LITAF","HCLS1","MT-ND4L","RPS26","CCL2","MT-ATP8","APOC1","MALAT1","DUSP1","HAVCR2","NUPR1","MT-ND4","EEF1A1","CCL7","SPP1","CXCL8","MT-ATP6","RPL39","FOS"},"Realization_2":{"S100A10","MT-CO2","SLC8A1","RPS26","RPL37A","D2HGDH","CCL2","MT-ATP8","BRI3","MALAT1","HLA-B","CXCL8","CCL7","NUPR1","MT-ND4","TSC22D3","APOC1","SRGN","CTSB","HLA-A"},"Realization_3":{"PSME2","RPS26","CCL2","MALAT1","MT-ATP8","APOC1","CCL7","MT-CO2","RPS28","CXCL8","MT-ATP6","SPP1","MT-ND4","NUPR1","HLA-DQA2","KLF6","HLA-A","APOBEC3A","RPL39","HLA-B"},"Realization_4":{"B2M","IRF7","HSBP1","UTRN","RPS26","CCL2","MALAT1","APOC1","MT-ATP8","CXCL8","S100A8","MT-ND4","TMSB4X","GAPDH","MT-ATP6","NUPR1","IGLV3-19","MX1","HLA-DQA2","FOS"},"Realization_5":{"EEF2","MT-CO2","RPS26","CCL2","MALAT1","CXCL8","HLA-B","DPAGT1","MRPS27","MT-ATP8","NUPR1","HLA-DQA2","MT-ND4","B2M","APOC1","TIMP1","RPL39","DUSP1","TXN","HLA-A"}}
#############################Plots#################################################################
font={'family':'normal','style':'normal','variant':'normal', 'weight':'normal','size':12}
plt.rc('font',**font)
# matplotlib.rc('font',**font)

def plot_confusion_matrix(cm, classes, ax, normalized=True, cmap=cmap): #bone
    ax = ax or plt.gca()
    norm_cm = cm
    if normalized:
        norm_cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    sns.set(font_scale=.9)
    sns.heatmap(norm_cm, annot=cm, fmt='g', xticklabels=classes, yticklabels=classes, cmap=cmap,cbar = False)


def shap_plot(shap_values, X,ax): #bone
    ax = ax or plt.gca()
    shap.summary_plot(shap_values, X)


fig = plt.figure(figsize=(55,55))
grid = plt.GridSpec(6,7,wspace=0.4,hspace=0.3)
# Confusion matrix test
ax = plt.subplot(grid[0:2,0:2])
cm_display = ConfusionMatrixDisplay(confusion_matrix=cm,display_labels=['Moderate','Severe']).plot(cmap=cmap, ax=ax,colorbar=False)
# Venn diagram
ax = plt.subplot(grid[0:3,2:7])
venn(lista,ax=ax,fontsize=8,legend_loc="upper right")
# shap values
ax = plt.subplot(grid[3:6,4:7])
ax.tick_params(labelsize=1)
shap_plot(shap_values, X,ax)

# confusion matrix validation cm_val
ax = plt.subplot(grid[3:6,0:2])
cm_display = ConfusionMatrixDisplay(confusion_matrix=cm_val,display_labels=['Moderate','Severe']).plot(cmap=cmap, ax=ax,colorbar=False)


# Saving
plt.savefig('.../figures/figures2_paper_corrected.png',format='png', dpi=300, bbox_inches='tight')


##########################################################################################################
##########################################################################################################
##########################################################################################################
