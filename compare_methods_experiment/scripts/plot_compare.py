import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.metrics import accuracy_score, confusion_matrix


ddG = pd.read_csv('ddG_compare.csv')

def correlationplots(ddG):
    
    correlations = ddG.corr(method='pearson')
    corr = correlations.iloc[0][1:]
    
    y_list = list(corr.index)
    n = len(y_list)

    x = ddG['avg_ddG']
    hfont = {'fontname':'Arial'}
    
    fig, ax = plt.subplots(ncols=n, nrows=1, sharey=True, figsize=(18,5))
    fig.subplots_adjust(wspace=0, hspace=0)
    
    for y in range(len(y_list)):        
        ax[y].plot(x, ddG[y_list[y]], color='black', alpha=0.6, ls='None',marker='o',markersize=3.)    
        m, b = np.polyfit(x, ddG[y_list[y]], 1)
        ax[y].plot(x, corr[y_list[y]]*x+b, color = 'black')
        ax[y].set_title(f"r = {round(corr[y],3)}", **hfont)
        ax[y].set_xlabel('Experimental '+r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)', **hfont)
    
    fig.text(0.09, 0.5, 'Predicted '+r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)', **hfont, va='center', rotation='vertical')
    #fig.set_ylabel('Predicted '+r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)', **hfont)
    #plt.show()
    plt.savefig("regplot.pdf", bbox_inches='tight')

    return

correlationplots(ddG)

#/data/user/shared_projects/p53_jmb_2021/compare_methods_experiment

def classification_columns(ddG):
    
    for i in range(1,len(ddG.columns)):
        col = ddG.columns[i]
        conditions = [ ddG[col] <= -1.2, (ddG[col] >= -1.2) & (ddG[col] <= 1.2), (ddG[col] > 1.2) & (ddG[col] <= 3), ddG[col] > 3 ]
        choices = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
        ddG[ddG.columns[i]+"_classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG

classification_columns(ddG)

def classification_plot(test_list, ddG):
    
    true = ddG['avg_ddG_classification']
    hfont = {'fontname':'Arial'}
    classnames = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
    fig, ax = plt.subplots(ncols=5, nrows=1, sharey=True, figsize=(18,4))
    fig.subplots_adjust(wspace=0, hspace=0)

    
    for i in range(len(test_list)-1):
        accuracy = accuracy_score(true, ddG[test_list[i]])   
        values = np.mat(confusion_matrix(true, ddG[test_list[i]]))
        df = pd.DataFrame(values, index=classnames, columns=classnames)        
        h = sns.heatmap(df, annot=True, cmap="viridis", cbar=False, ax=ax[i], vmax=21)
        ax[i].xaxis.set_ticklabels(h.xaxis.get_ticklabels(), rotation=90, ha='right')    
        #ax[i].set_xlabel('True label', **hfont)
        ax[i].set_title(f"Accuracy: {round(accuracy,3)}", **hfont)
        
    i = 4
    accuracy = accuracy_score(true, ddG[test_list[i]])   
    values = np.mat(confusion_matrix(true, ddG[test_list[i]]))
    df = pd.DataFrame(values, index=classnames, columns=classnames)        
    h = sns.heatmap(df, annot=True, cmap="viridis", vmax=21)
    ax[i].xaxis.set_ticklabels(h.xaxis.get_ticklabels(), rotation=90, ha='right')    
    #ax[i].set_xlabel('True label', **hfont)
    ax[i].set_title(f"Accuracy: {round(accuracy,3)}", **hfont)
        
    #fig.yaxis.set_ticklabels(ax[i].yaxis.get_ticklabels(), rotation=0, ha='right')
    fig.text(0.04, 0.5, 'Predicted label', **hfont, va='center', rotation='vertical')
    
    #plt.show()
    
    plt.savefig("heatmap.pdf", bbox_inches='tight')
    
    return    

test_list = list(ddG.columns[8:])
classification_plot(test_list, ddG)

