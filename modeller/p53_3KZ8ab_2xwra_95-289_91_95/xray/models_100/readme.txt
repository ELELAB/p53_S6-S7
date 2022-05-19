Read me for preparation of 3KZ8 dimer from tetramer complex analysis. 
=======================================================================

The input file 3KZ8 was was fetched into a pymol session where:
        1) Water and iodoine was removed while ZN atoms were kept
        2) Symmetry mates for the protein and DNA were created
        3) Chain A and B of the protein was exported as a molecule in 3kz8_01.pdb 

Part 1: Pymol
=====================
3kz8_01.pdb and the template 2XWR_clean_A.pdb were opened into a pymol session where:
        1) The reference structure 2XWR_clean_A.pdb was aligned to the Ca's of 3kz8_01.pdb chain A:
           select calpha3KZ8_A,  (name CA) and  chain A and ss S and 3kz8_01
           align 2XWR_clean_A, calpha3KZ8_A
           save the aligned 2XWR_clean_A as a molecule called 2XWRa_aligned_A.pdb
        2) The reference structure 2XWR_clean_A.pdb was aligned to the Ca's of 3kz8_01.pdb chain B:
           select calpha3KZ8_B,  (name CA) and  chain B and ss S and 3kz8_01
           align 2XWR_clean_A, calpha3KZ8_B
           save the aligned 2XWR_clean_A as a molecule called 2XWRa_aligned_B.pdb

        
Part 2: Renaming and Upload
============================
The Zinc atoms in 3kz8_01.pdb were renamed following 2XWR convention, residue number 300 and moved after TER entry of the corresponding chain
 
Upload to the server the pse session and 2XWRa_aligned_A.pdb and 2XWRa_aligned_B.pdb
(In this case a copy of 2XWRa_aligned.pdb called 2XWRa_aligned_2.pdb was created because two areas need modeling using the 2XWR template)

Part 2: Update the scripts
============================
Include in the py file of modeller modeller-multiple3.py the correct range of missing residues and missing heavy atoms in your pdbs
update the correct atoms for the restraints
update the xxx.ali file with the right name of the template pdbs 


Part 4: Running modeller to fix described issues and build missing items from template
====================================================================================
This requires a number of files: 

        1) xxx.ali 
           An alignment file for the areas that needs to be modelled. 
        2) modeller-multiple3.py
           A script that runs the modeller where: 
                2.1 modelling of tail: ('1:A', '6:A')
                2.2 Modelling of missing region: ('95:A','97:A') 
                2.3 modelling of missing heavy atoms: self.atoms['NZ:11:A']...
                2.4 restraints for modelling. 
           Running: tsp -N 1 mod10.1 modeller-multiple3.py

Part 5: Clean up and prepare for analysis
====================================================================================
mkdir models
mv *B9999* models
rm *V9999*
rm *D000* 
