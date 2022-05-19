import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("/data/user/shared_projects/p53_jmb_2021/foldX_DNA/results/ddG_overview.csv")

ddG_0 = df[df.Binding_energy == 0]
print(f"Number of mutations without any change is {len(ddG_0)}")
#656 mutations does not have any effect on DNA binding

ddG_c = df[df.Binding_energy != 0]
print(f"Number of mutations with any change is {len(ddG_c)}")

def classification(ddG):
    
    conditions = [ ddG.Binding_energy <= -1.2, (ddG.Binding_energy >= -1.2) & (ddG.Binding_energy <= 1.2), (ddG.Binding_energy > 1.2) & (ddG.Binding_energy <= 3), ddG.Binding_energy > 3 ]
    choices = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
    ddG["classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG

ddG_c = classification(ddG_c)

ddG_effect = ddG_c[ddG_c.classification != 'neutral']

def plot(ddG, w, figname):
    
    hfont = {'fontname':'Arial'}
    
    ddG = ddG.sort_values(by='Mutation', key=lambda col: col.str[2:-1].astype(int))

    #x and y values
    x = ddG['Mutation']
    y = ddG['Binding_energy']
    
    fig, ax = plt.subplots(figsize=(w,5))
    fig.subplots_adjust(wspace=0, hspace=0)
    
    #plotting
    ax.scatter(x, y, c='black', alpha=0.6, ls='None', marker='o') 
    ax.set_xlabel('Mutations', **hfont)    
    ax.set_xticklabels(x, rotation=90, ha='right')
    ax.margins(x=0.01)
        
    ax.set_ylabel(r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)', **hfont)
        
    plt.savefig(f"scatterplot_DNA_{figname}.pdf", bbox_inches='tight')

    return

plot(ddG_c, 40, "non_0")

plot(df, 100, "all")

plot(ddG_effect, 10, "effect")
