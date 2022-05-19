import pandas as pd
import numpy as np

path = "."
file = "saturation_mutagenesis_overview.csv"
file_cancermuts = 'cancermuts_mutations.txt'

#def cancer_muts_import(file):
    
#    df = pd.read_csv(file) 
#    mutation_list = []

#    wt = list(df['WT.residue'])
#    pos = list(df['prot.Position'])
#    mut = list(df['Mutated.residue'])
    
#    for i in range(len(df)):
#        if pos[i] > 93 and pos[i] < 291:
#            mutation = wt[i]+str(pos[i])+mut[i]
#            mutation_list.append(mutation)
#    
#    return mutation_list

def read_cancermuts(file):       
    mutation_list = [] 
    cancermuts = open(file, "r")
    content = cancermuts.readlines()
    for line in content:
        mutation_list.append(line[:-2])
    return mutation_list

def calculate_average(df):

    mutations = []
    standard_deviation = []
    average = []
    
    for i in range(len(df)):
        
        mutations.append(df.mutant[i])
    
        std = np.std([df.ddG_rosetta_exp[i], df.ddG_mutatex_exp[i], df.ddG_mutatex_pdb[i], df.ddG_mutatex_md[i]])
        standard_deviation.append(std)
        
        #avg = np.average([df.ddG_rosetta_exp[i], df.ddG_mutatex_exp[i], df.ddG_mutatex_pdb[i], df.ddG_mutatex_md[i]])
        #average.append(avg)
        
        if std < 1.2:
            avg = np.average([df.ddG_rosetta_exp[i], df.ddG_mutatex_exp[i], df.ddG_mutatex_pdb[i], df.ddG_mutatex_md[i]])
            average.append(avg)
            
        else:
            avg = df.ddG_mutatex_md[i]
            average.append(avg)
    
    d = {'mutant': mutations, 'avg': average, 'std': standard_deviation}
    
    df = pd.DataFrame(data=d)
    
    return df

def classification_columns(ddG):
    
    conditions = [ ddG['avg'] <= -1.2, (ddG['avg'] >= -1.2) & (ddG['avg'] <= 1.2), (ddG['avg'] > 1.2) & (ddG['avg'] <= 3), ddG['avg'] > 3 ]
    choices = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
    ddG["classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG   

def extract_loops(df):

    placement = []

    for i in df.mutant:
        if int(i[1:-1]) in range(114,124):
            placement.append("L1")
        elif int(i[1:-1]) in range(207,213):
            placement.append("S6-7")
        else:
            placement.append("N")
            
    df['loop'] = placement
    
    L1_df = df[df.loop == 'L1']
    S67_df = df[df.loop == 'S6-7']
    
    return L1_df, S67_df

def extract_cancermuts(df, cancermuts_list):
    
    placement = []
    
    for i in df.mutant:
        if i in cancermuts_list:
            placement.append("cancermuts")
        else:
            placement.append("N")
    
    df['cancermuts'] = placement
    
    cancermuts_df = df[df.cancermuts == 'cancermuts']
    
    cancermuts_df = cancermuts_df.drop(columns=['cancermuts'])
    
    return cancermuts_df
        
    
df_initial = pd.read_csv(path+file, index_col=0)
avg_df = calculate_average(df_initial)    
df = classification_columns(avg_df)  
#l1_df, s67_df = extract_loops(df)
cancermuts = read_cancermuts(path+file_cancermuts)
#cancermuts = cancer_muts_import(file_cancermuts)
cancermuts_df = extract_cancermuts(df, cancermuts)


#l1_df.to_csv('L1_mutations.csv')
#s67_df.to_csv('S6-S7_mutations.csv')
cancermuts_df.to_csv('cancermuts_mutations.csv')
avg_df.to_csv('saturation_mutagenesis_avg.csv')
