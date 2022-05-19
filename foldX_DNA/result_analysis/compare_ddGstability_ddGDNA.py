#compare_ddG_DNA, DDG stab
import pandas as pd
import numpy as np

def DNA_classification(ddG):
    
    conditions = [ ddG.Binding_energy <= -0.8, (ddG.Binding_energy >= -0.8) & (ddG.Binding_energy <= 0.8), ddG.Binding_energy > 0.8]
    choices = [ "stabilizing", 'neutral', 'destabilizing']
    
    ddG["DNA_classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG

def stability_classification(ddG):
    
    conditions = [ ddG.Average <= -1.2, (ddG.Average >= -1.2) & (ddG.Average <= 1.2), (ddG.Average > 1.2) & (ddG.Average <= 3), ddG.Average > 3 ]
    choices = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
    ddG["stability_classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG

#DNA ddG
def import_transform_DNA(file):
    
    #import csv
    df_DNA = pd.read_csv(file)
    
    #rename df_DNA to exclude chain name 
    mutant = []
    for i in range(len(df_DNA)):
        new_name = df_DNA.Mutation[i][0] + df_DNA.Mutation[i][2:]
        mutant.append(new_name)
    
    #add new names
    df_DNA['mutant'] = mutant
    df_DNA = DNA_classification(df_DNA)
    #Remove all netral interactions
    #df_DNA = df_DNA[df_DNA.classification != "neutral"]
    
    return df_DNA


#import the selected neutral mutations
def import_neutral_mutations(file):
    
    sel_neu = open(file, "r")
    sel_neu = sel_neu.read()
    sel_neu = sel_neu.split("\n")
    sel_neu = sel_neu[1:-1] #remove header and last newline.
    
    for i in range(len(sel_neu)):
        if sel_neu[i][-1] == "*":
            sel_neu[i] = sel_neu[i][0:-1] #269
    
    return sel_neu

def import_stability_calc(file):

    df_all = pd.read_csv(file, index_col=0)
    df_all = df_all.drop(['class_rosetta_exp', 'ddG_rosetta_pdb','class_rosetta_pdb', 'class_mutatex_exp','class_mutatex_pdb','class_mutatex_md'], axis=1)
    
    mean = []
    deviation = []
    
    for i in range(len(df_all)):
        if abs(df_all.loc[i][1:].mean()-df_all.loc[i].ddG_mutatex_md) < 1.2:
            value = df_all.loc[i][1:].mean() 
            std = df_all.loc[i][1:].std()
            mean.append(value)
            deviation.append(std)
        else:
            value = df_all.loc[i].ddG_mutatex_md
            std = df_all.loc[i][1:].std()
            mean.append(value)
            deviation.append(std)
    
    df_all['Average'] = mean
    df_all['std'] = deviation
    
    df_all = stability_classification(df_all)
    
    return df_all


def build_analysis_df(df_all, df_DNA, selected_mutations):
    #create a dataframe with the selected neutral mutations and their DNA ddG value
    ddG_all_sel_neu = df_all[df_all.mutant.str.contains('|'.join(selected_mutations))]
    ddG_DNA_sel_neu = df_DNA[df_DNA.mutant.str.contains('|'.join(selected_mutations))]
    
    #combine on "mutant"
    analysis_df = ddG_all_sel_neu.merge(ddG_DNA_sel_neu, on='mutant', how='inner')
    #clean up df
    analysis_df = analysis_df.drop(['Mutation'], axis=1)
    #exclude neutral DNA ddGs

    analysis_df = analysis_df[analysis_df.DNA_classification != "neutral"]
    
    return analysis_df


selected_mutations = import_neutral_mutations("selected_neutral_list.txt")
df_DNA = import_transform_DNA("ddG_overview.csv")
df_all = import_stability_calc("saturation_mutagenesis_overview.csv")

neu_ddG_DNA_disrupt = build_analysis_df(df_all, df_DNA, selected_mutations)

export_df = neu_ddG_DNA_disrupt.drop(['ddG_rosetta_exp', 'ddG_mutatex_exp', 'ddG_mutatex_pdb','ddG_mutatex_md'], axis=1)

export_df.to_csv('table_s2.csv')
