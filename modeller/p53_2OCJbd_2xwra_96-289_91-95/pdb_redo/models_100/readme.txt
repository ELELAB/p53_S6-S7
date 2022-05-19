#As template we here used the X-ray structure of p53 DBD homodimer 2OCJ using the orientation defined by chain B and D 
#We here used the final structure optimized from pdb_redo
#see repository in pdbs_data for further info

cp ../../../../../../../pdbs_data/tcga_3d/p53/zinc_bound/pdbs_redo_processed/split_chain/2ocj_final_tot_clean_BD.pdb .

#We here included p53 DBD residues 96-289 from the PDB entry 2OCJ
#We used MODELLER to reconstruct only the missing region 91-95 using as template the structure of p53 DBD PDB 2XWR chain A 
cp ../../../../../../../pdbs_data/tcga_3d/p53/zinc_bound/pdbs_processed/split_chain/2XWR_clean_A.pdb .

#while all the rest of the structure is kept frozen (no optimization by MODELLER) during the modeling
#We used pymol to align the template 2XWR_clean_A.pdb to chain B and chain D of 2OCJ and saved the coordinates as 2XWR_clean_A_alignB.pdb and 2XWR_clean_A_alignD.pdb. For the alignment we used all the C alpha atoms of the secondary structures in p53 DBD   
#We converted the template structure to  generate the alignment file (called xxx.ali) 
#We modified the alignment file (xxx.ali) to include the templates
#We modified the MODELLER run script named modeller-multiple3.py to include the correct selection of the atoms/residues that shouldn’t or should be optimized during the modeling
#We reconstructed and optimized residues 91-96 using six harmonic restraints to preserve the distance between atoms known to form interactions in 2XWR described in Natan E et al.. J. Mol Biol. 2011 Jun 10;409(3):358-68 PMID:21457718
#We used mod10.1 to predict the 100 models with the command 

mod10.1 modeller-multiple3.py

#we included all the models in the folder models

mkdir models

#we moved in the models folder the final models and removed the partial files created by MODELLER

mv *B9999* models
rm *V999*
rm *D000*

#We performed the selection of models in the modeller_analysis folder

/data/user/shared_projects/p53_jmb_2021/modeller_analysis/zinc-bound/p53_2OCJbd_2xwra_96-289_91-95/pdb_redo/models_100
