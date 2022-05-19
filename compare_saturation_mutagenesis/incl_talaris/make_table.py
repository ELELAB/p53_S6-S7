#rosetta

import pandas as pd
import numpy as np
import glob

def classification_columns(ddG):
    
    conditions = [ ddG['total_score'] <= -1.2, (ddG['total_score'] >= -1.2) & (ddG['total_score'] <= 1.2), (ddG['total_score'] > 1.2) & (ddG['total_score'] <= 3), ddG['total_score'] > 3 ]
    choices = [ "stabilizing", 'neutral', 'destabilizing', 'highly destabilizing' ]
    
    ddG["classification"] = np.select(conditions, choices, default=np.nan)
    
    return ddG

def pick_up_rosetta(path, name):
    
    df = pd.read_csv(path)
    df = df.loc[df['state'] == 'ddg']
    df = df[['mutation_label', 'total_score']]
    df = classification_columns(df)
    
    df = df.rename(columns={"mutation_label": "mutant", "total_score": f"ddG_{name}", "classification": f"class_{name}"})
    
    return df
    
def pick_up_mutatex(path, name):
    mutatex_aa_order = ['G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H'] 
    files = glob.glob(f'{path}/*')
    
    substitution_list = []
    value_list = []
    
    for file in files:
        filedf = pd.read_csv(file)
        values = filedf['# avg\tstd\tmin\tmax'] 
        for i in range(len(values)):
            t = file.split("/")[-1]
            t = t[0]+t[2:]
            AA_sub = f"{t}{mutatex_aa_order[i]}"
            substitution_list.append(AA_sub)
            avg = float(values[i].split(" ")[0])
            value_list.append(avg)
            
    d = {'mutant': substitution_list, 'total_score': value_list}
    
    df = pd.DataFrame(data=d)
    
    df = classification_columns(df)
    
    df = df.rename(columns={"total_score": f"ddG_{name}", "classification": f"class_{name}"})
    
    return df
    

p = "../"

r_exp = pick_up_rosetta(p+"rosetta_analysis/zinc_bound/2XWRa_91-289/exp/ref2015/aggregate/ddg_mutations_aggregate.csv", "rosetta_exp")

r_pdb_redo = pick_up_rosetta(p+"rosetta_analysis/zinc_bound/2XWRa_91-289/pdbredo/ref2015/aggregate_march2022/ddg_mutations_aggregate.csv", "rosetta_pdb")

all_ddG = pd.merge(r_exp, r_pdb_redo, on='mutant', how='left')

r_talaris_exp = pick_up_rosetta("/data/raw_data/computational_data/rosetta_data/tcga_3d/p53/zinc_bound/2XWRa_91-289/exp/talaris2014/output_aggregate/ddg_mutations_aggregate.csv", "talaris_exp")

all_ddG = pd.merge(all_ddG, r_talaris_exp, on='mutant', how='left')

m_exp = pick_up_mutatex(p+"mutatex_analysis/zinc_bound/2XWRa_91-289/exp_noHOH/saturation_new/final_averages", "mutatex_exp")
all_ddG = pd.merge(all_ddG, m_exp, on='mutant', how='left')

m_pdb = pick_up_mutatex(p+"mutatex_analysis/zinc_bound/2XWRa_91-289/pdbredo_noHOH/saturation/final_averages/", "mutatex_pdb")
all_ddG = pd.merge(all_ddG, m_pdb, on='mutant', how='left')

m_md = pick_up_mutatex(p+"mutatex_analysis/zinc_bound/2XWRa_91-289/md/replicate1/CHARM22star/saturation_rename_NEW/final_averages/", "mutatex_md")
all_ddG = pd.merge(all_ddG, m_md, on='mutant', how='left')

all_ddG.to_csv('saturation_mutagenesis_overview.csv')
