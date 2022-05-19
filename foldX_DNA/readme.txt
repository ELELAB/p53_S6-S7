Readme 
=============================================================

#To run the snakemake pipeline, python 3.7 or above needs to be activated. 

module load python/3.7/modulefile 

#The pipeline needs two inputs, a PDB structure and a list of mutations. 

#PDB structure:
#==================
#The PDB structure is pdb ID: 3KZ8. A p53 structure that is in complex with DNA, which has undergone modelling with modeller based on template PDB #ID: 2XWR. 

#List of mutations:
#==================
    
#   from /data/user/shared_projects/p53_jmb_2021/cancermuts/jan2022/, metatable.csv

#used to create the individuals list using create_individual_list.py. Notice that the script is updated based 
#on the positions available in the PDB structure. 

#Input data:
#	create_individual_list.py
#	metatable.csv

cd input_data
python3 create_individual_list.py

#Input data:
#	create_individual_list.py
#	individual_list.txt
#	metatable.csv

cd ..
cp input_data/individual_list.txt .

# config.yaml:
#==================

#update as needed


# snakefile:
#==================

#snakemake pipeline was run on individuals list by activating the python 3 module with:

tsp -N 4 snakemake --cores 4

