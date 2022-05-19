#Compare mutatex, rosetta 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter

path = "."
file = "saturation_mutagenesis_overview.csv"

#describe the agreement

def kde_plot(path, file):
    
    df = pd.read_csv(path+file, index_col=0)

    sns.set(style="white")
    
    fig, axs = plt.subplots(5, 1, figsize=(5, 20))
    
    s = sns.kdeplot(data=df, x="ddG_rosetta_exp", color="black", alpha=0.7, ax=axs[0])
    s.set_xlim(-7, 30)
    s.set_ylim(0, 0.32)
    s.set_xlabel('Rosetta, X-ray structure')
    
    d = sns.kdeplot(data=df, x="ddG_rosetta_pdb", color="black", alpha=0.7, ax=axs[1])
    d.set_xlim(-7, 30)
    d.set_ylim(0, 0.32)
    d.set_xlabel('Rosetta, PDBredo structure')
    
    f = sns.kdeplot(data=df, x="ddG_mutatex_exp", color="black", alpha=0.7, ax=axs[2])
    f.set_xlim(-7, 30)
    f.set_ylim(0, 0.32)
    f.set_xlabel('MutateX, X-ray structure')
    
    g = sns.kdeplot(data=df, x="ddG_mutatex_pdb", color="black", alpha=0.7, ax=axs[3])
    g.set_xlim(-7, 30)
    g.set_ylim(0, 0.32)
    g.set_xlabel('MutateX, PDBredo structure')
    
    h = sns.kdeplot(data=df, x="ddG_mutatex_md", color="black", alpha=0.7, ax=axs[4])
    h.set_xlim(-7, 30)
    h.set_ylim(0, 0.32)
    h.set_xlabel('MutateX, Molecular Dynamics Ensemble')
    
    plt.savefig("KDE.pdf", bbox_inches='tight')
    
    plt.show()
    
    return fig

def find_majority(path, file):
    
    df = pd.read_csv(path+file, index_col=0)
    
    classification_majority = []
    classification_fraction = []
    
    for i in range(len(df)):
        line  = df.loc[i]
        classes = [line.class_rosetta_exp, line.class_rosetta_pdb, 
                   line.class_mutatex_exp, line.class_mutatex_pdb, 
                   line.class_mutatex_md]
         
        vote_count = Counter(classes)
        top_five = vote_count.most_common(5)
        
        classification_majority.append(top_five[0][0])
        classification_fraction.append(top_five[0][1]/5)
    
    df['class_majority'] = classification_majority
    df['class_fraction'] = classification_fraction
    
    return df

def find_percentage_agreement(classification_df):
    
    classes = ['destabilizing', 'highly destabilizing', 'neutral', 'stabilizing']
    fractions = [0.4, 0.6, 0.8, 1.0]
    
    classification =[]  
    value =[]
    counting = []
    
    
    for i in classes:        
        for j in fractions:
            
            count = len(classification_df[(classification_df.class_majority==i) & (classification_df.class_fraction==j)])/len(classification_df[(classification_df.class_majority==i)])
            
            classification.append(i)  
            value.append(j)
            counting.append(count)
        
   
    d = {"classifi": classification, 'method agreement, fraction': value, 'volumne of agreement, fraction': counting}
    
    df = pd.DataFrame(d)
    
        
    return df

def fraction_plot(path, file):
    
    classification_df = find_majority(path, file)
    percentage_df = find_percentage_agreement(classification_df)
    
    g = sns.lineplot(data=percentage_df, x="method agreement, fraction", y="volumne of agreement, fraction", hue="classifi", markers=True, dashes=False)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4)
    
    plt.savefig("fractions.pdf", bbox_inches='tight')
    plt.show()
    
    return g

kde_plot(path, file)
fraction_plot(path, file)

