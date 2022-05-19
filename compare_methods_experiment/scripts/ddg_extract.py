import pandas as pd

#MutateX functions

def doc_names(exp, chain):
    
    #extract all mutation
    Mut_list = list(exp.Mutation)
    
    #create an empty list of files where these can be found
    file_list = []
    AA = []
    
    #looping over each mutation and annotating the file
    for mutation in range(len(Mut_list)):
        
        mutation_file = Mut_list[mutation][0] + chain + Mut_list[mutation][1:-1]
        file_list.append(mutation_file)
        
        mut = Mut_list[mutation][-1]
        AA.append(mut)
    
    return file_list, AA


def avg_ddg_mutatex(relevant_documents):
    mutatex_aa_order = ['G','A','V','L','I','M','F','W','P','S','T','C','Y','N','Q','D','E','K','R','H']    
    avg = [] 
    #std = []
    for i in range(len(relevant_documents[0])):
        file = pd.read_csv(relevant_documents[0][i]) 
        values = file['# avg\tstd\tmin\tmax'] 
        line = mutatex_aa_order.index(relevant_documents[1][i]) 
        value = values[line].split(" ") 
        avg_value = float(value[0]) 
        #std_value = float(value[1])
        avg.append(avg_value)
        #std.append(std_value)
    return avg #, std

#rosetta Function

def avg_ddg_rosetta(mutation_list, table):
    #rosetta does not offer STD values readily, to ask Vale
    file = pd.read_csv(table) 
    avg = [] 
    #std = []
    for i in range(len(mutation_list)):
        file[file['mutation_label'] == mutation_list[i]] 
        line = file.loc[(file['mutation_label'] == mutation_list[i]) & (file['state']=='ddg')] 
        ddg = float(line.total_score) 
        avg.append(ddg)
        #deviation = float(line.??)
        #std.append(deviation)
    return avg #, std
