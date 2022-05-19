from Bio.PDB.Polypeptide import one_to_three
import pandas as pd

mutation_files = [ 'germline_from_literature.txt',      'germline_p53DB_list.txt',      'gnomAD_low_freq.txt',      'somatic_p53DB_list.txt' ]
out_files      = [ 'germline_from_literature_HGVS.csv', 'germline_p53DB_list_HGVS.csv', 'gnomAD_low_freq_HGVS.csv', 'somatic_p53DB_list_HGVS.csv' ]


for fi, fo in zip(mutation_files, out_files):
    with open(fi, 'r') as fhi:
        mutations = []
        for line in fhi.read().splitlines():
            res1 = one_to_three(line[0])
            res1 = res1[0] + str.lower(res1[1:])
            res2 = one_to_three(line[-1])
            res2 = res2[0] + str.lower(res2[1:])
            mutations.append(f"p.{res1}{line[1:-1]}{res2}")
        type_col = ['mutation'] * len(mutations)
        function_col = [''] * len(mutations)
        df = pd.DataFrame({ 'name' : type_col,
                            'type' : type_col,
                            'site' : mutations,
                            'function' : function_col,
                            'reference' : function_col })
        df.to_csv(fo, sep=';', index=False)
