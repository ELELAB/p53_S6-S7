import pandas as pd

df = pd.read_csv('metatable.csv')

wt = list(df['WT residue'])
pos = list(df['Position'])
mut = list(df['Mutated residue'])

with open('cancermuts_mutations.txt', 'w') as mut_file:
    for i in range(len(df)):
        if pos[i] > 90 and pos[i] < 290 and mut[i] != "nan":
            mut_file.write(wt[i]+str(pos[i])+mut[i]+'\n')
