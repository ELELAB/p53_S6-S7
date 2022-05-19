import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

L1 = pd.read_csv('L1_mutations.csv', index_col=0)
s67 = pd.read_csv('S6-S7_mutations.csv', index_col=0)
cancermuts = pd.read_csv('cancermuts_mutations.csv', index_col=0)

def plot_loops(loop_df, name):
    
    x = loop_df['mutant']
    y = loop_df['avg']
    
    plt.figure(figsize=(15,4))
    plt.errorbar(x, y, yerr=loop_df['std'], fmt="o", color='black', alpha=0.3)
    
    categories = np.unique(loop_df['classification'])
    colors = np.linspace(0, 1, len(categories))
    colordict = dict(zip(categories, colors)) 
    loop_df["Color"] = loop_df['classification'].apply(lambda x: colordict[x])
    
    plt.scatter(x, y, c=loop_df["Color"], alpha=0.9)
    
    plt.xticks(rotation=90)
    
    plt.xlabel('Mutations')
    plt.ylabel('Predicted ΔΔG (kcal/mol)')
    
    plt.savefig(f"{name}_overview.pdf", bbox_inches="tight")
    
    plt.show()
    
    return

plot_loops(L1, 'L1')
plot_loops(s67, 'S6-S7')

def plot_cancermuts(df):
    
    x = df['mutant']
    y = df['avg']
    
    plt.figure(figsize=(50,4))
    plt.errorbar(x, y, yerr=df['std'], fmt="o", color='black', alpha=0.3)
    
    categories = np.unique(df['classification'])
    colors = np.linspace(0, 1, len(categories))
    colordict = dict(zip(categories, colors)) 
    df["Color"] = df['classification'].apply(lambda x: colordict[x])
    
    plt.scatter(x, y, c=df["Color"], alpha=0.9)
    
    plt.xticks(rotation=90)
    
    plt.xlabel('Mutations')
    plt.ylabel('Predicted ΔΔG (kcal/mol)')
    
    plt.savefig("cancermuts_overview.pdf", bbox_inches="tight")
    
    plt.show()
    
    return

plot_cancermuts(cancermuts)

