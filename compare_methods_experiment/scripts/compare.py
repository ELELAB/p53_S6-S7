import pandas as pd
import os
from ddg_extract import doc_names, avg_ddg_mutatex, avg_ddg_rosetta

###############################################################################
#   Prepare the path and file
###############################################################################

#define path
path = os.getcwd()

#import comparison file
exp = pd.read_csv("tableExport_DBD_single.csv")

#Extract relevant columns
exp = exp[['Mutation', 'ddG']]

#reverse ddG value so it is compatible with mutateX results.
exp['ddG'] = exp['ddG']*-1

avg = exp.groupby('Mutation').mean()['ddG']
#std = exp.groupby('Mutation').std()['ddG']

exp = avg.to_frame()

exp = exp.rename(columns={'ddG' : 'avg_ddG'})

exp = exp.fillna(0)

exp = exp.reset_index()

###############################################################################
#   Running the collection of data
###############################################################################


#   MutateX
###############################################################################
#a list used to find the files created by mutateX
documents = doc_names(exp, "A")

os.chdir(path+"/md_mutatex/final_averages/")
avg1 = avg_ddg_mutatex(documents)
exp["ddG_md_mutateX"] = avg1
#exp["std_md_mutateX"] = avg1[1]

os.chdir(path+"/xray_mutatex/final_averages/")
avg2 = avg_ddg_mutatex(documents)
exp["ddG_xray_mutateX"] = avg2
#exp["std_xray_mutateX"] = avg2[1]

os.chdir(path+"/pdbredo_mutatex/final_averages/")
avg3 = avg_ddg_mutatex(documents)
exp["ddG_pdbredo_mutateX"] = avg3
#exp["std_pdbredo_mutateX"] = avg3[1]

#   Rosetta
###############################################################################

os.chdir(path+"/xray_rosetta/")
avg4 = avg_ddg_rosetta(list(exp['Mutation']), "ddg_mutations_aggregate.csv")
exp["xray_rosetta"] = avg4

os.chdir(path+"/pdbredo_rosetta/")
avg5 = avg_ddg_rosetta(list(exp['Mutation']), "ddg_mutations_aggregate.csv")
exp["pdbredo_rosetta"] = avg5


#   Export the file
###############################################################################
os.chdir(path)
exp.to_csv("ddG_compare.csv", index=False)




