import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# color for our classes of consequence on protein stability
colors = {'stab'    : "#440154",
          'neut'    : "#31688e",
          'lowdes'  : "#35b779",
          'highdes' : "#fde725"}

# function to color according to stability
def ddg_to_color(x, col):
    if x[col] < -1.2:
        return colors['stab']
    if -1.2 <= x[col] <= 1.2:
        return colors['neut']
    if 1.2 < x[col] <= 3:
        return colors['lowdes']
    if x[col] > 3:
        return colors['highdes']

# parse data
ddgs = pd.read_csv('../../saturation_mutagenesis_overview.csv').set_index('mutant')
mutations = pd.read_csv('../../cancermuts_mutations.txt', header=None)[0]

# filter dataframe keeping only cancer mutations
selected_ddgs = ddgs[ ddgs.index.isin(mutations)]

# process to create columns we need
selected_ddgs['positions'] = selected_ddgs.apply(lambda x: x.name[:-1], axis=1)

# split in mutatex and rosetta datasets
data_mutatex = selected_ddgs.groupby('positions').ddG_mutatex_md.agg(list).explode()
data_rosetta = selected_ddgs.groupby('positions').ddG_rosetta_pdb.agg(list).explode()

# pair sequence position (number) with number+residue
sorted_positions = sorted(data_mutatex.index.unique(), key=lambda x: int(x[1:]))
xs = dict( [ (s, int(s[1:])) for s in sorted_positions ] )
#xs = dict(zip(sorted_positions, np.arange(1, len(sorted_positions)+1)))

# replace index with just residue number
data_mutatex = data_mutatex.rename(index=xs).to_frame()

# add color column according to classification
data_mutatex['color'] = data_mutatex.apply(ddg_to_color, col='ddG_mutatex_md', axis=1)

# to numpy array for plotting
data_mutatex_ar = data_mutatex['ddG_mutatex_md'].reset_index().to_numpy().T

# repeat for rosetta dataset
data_rosetta = data_rosetta.rename(index=xs).to_frame()
data_rosetta['color'] = data_rosetta.apply(ddg_to_color, col='ddG_rosetta_pdb', axis=1)
data_rosetta_ar = data_rosetta['ddG_rosetta_pdb'].reset_index().to_numpy().T

# plotting. Create subplots
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(7.48, 5))
fig.subplots_adjust(hspace=0)

# change font to Arial
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'

# customize plot
ax1.set_xticks([91] + list(np.arange(100, 290, 10)) + [289])
ax1.set_yticks(np.arange(-5, 30, 5))
ax2.set_yticks(np.arange(-5, 30, 5))
ax1.tick_params(labelsize=8)
ax2.tick_params(labelsize=8)

ax2.set_xlabel('p53 DBD sequence position')

ax1.set_xlim((90.5, 289.5))
ax1.set_ylim((-5, 30))
ax2.set_ylim((-5, 30))

# plotting
ax1.set_ylabel(r'FoldX $\Delta{}\Delta{}G$ (kcal/mol)')
ax1.scatter(data_mutatex_ar[0], data_mutatex_ar[1],
            c=data_mutatex['color'],
            s=3.5)

ax2.set_ylabel(r'Rosetta $\Delta{}\Delta{}G$ (kcal/mol)')
ax2.scatter(data_rosetta_ar[0], data_rosetta_ar[1],
		    c=data_rosetta['color'],
            s=3.5)

# add threshold line
ax1.axhline(y=1.2, xmin=-1, xmax=300, color='red', linestyle='--', lw=1)
ax2.axhline(y=1.2, xmin=-1, xmax=300, color='red', linestyle='--', lw=1)

# save figure
fig.savefig('figure_4b.pdf')

