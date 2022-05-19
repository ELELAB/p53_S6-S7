import pandas as pd
import matplotlib.pyplot as plt

#import comparison file
exp = pd.read_csv("tableExport_DBD_single.csv")
#Extract relevant columns
exp = exp[['Mutation', 'ddG']]
#reverse ddG value so it is compatible with mutateX results.
exp['ddG'] = exp['ddG']*-1

def plot_average(exp):
    
    avg = exp.groupby('Mutation').mean()['ddG']
    exp1 = avg.to_frame()
    exp1 = exp1.rename(columns={'ddG' : 'avg_ddG'})
    exp1 = exp1.fillna(0)
    exp1 = exp1.reset_index()
    
    plt.figure(figsize=(13,5))
    
    plt.plot(exp1['Mutation'], exp1["avg_ddG"], color='darkblue', alpha=0.6, ls='None',marker='o',markersize=3.)    
    plt.xticks(rotation=90)
    plt.ylim([-2, 5])
    plt.xlabel('Mutations')
    plt.ylabel('Average Experimental '+r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)')
    
    plt.savefig("avg_themomut_ddg.pdf", bbox_inches='tight')
    #plt.show()
    
    return

def plot_difference(exp):
    
    t = pd.concat(g for _, g in exp.groupby("Mutation") if len(g) > 1)

    plt.figure(figsize=(5,5))
    
    plt.plot(t['Mutation'], t["ddG"], color='darkblue', alpha=0.6, ls='None',marker='o',markersize=3.)    
    plt.xticks(rotation=90)
    plt.ylim([-2, 5])
    plt.xlabel('Mutations')
    plt.ylabel('Experimental '+r'$\Delta$'+r'$\Delta$'+'G (kcal/mol)')
    
    plt.savefig("differencies_in_thermomuts.pdf", bbox_inches='tight')
    #plt.show()
    
    return

plot_average(exp)
plot_difference(exp)
