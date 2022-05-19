Read me for preparation of 3KZ8 chain B for DNA complex analysis. 
=======================================================================

Part 1: Pymol
=====================
3KZ8 was fetched into a pymol session where:
        1) Water and iodoine was removed. 
        2) Symmetry mates for the DNA was created 
        3) Chain B of the protein and the DNA was exported as a molecule. 
        4) The reference structure 2XWR chain A was aligned to the Ca's:
                select calpha3KZ8,  (name CA) and 3KZ8_DNA_BCD and ss S
        
Part 2: Renaming and Upload
============================
3KZ8_DNA_BCD.pdb was opened as txt document and the symmetry mate was renamed as chain D. Zinc was renamed following 2XWR convention, number 300.
 
Upload of 3KZ8_DNA_BCD.pdb and 2XWRa_aligned.pdb
(In this case a copy of 2XWRa_aligned.pdb called 2XWRa_aligned_2.pdb was created because two areas need modeling using the 2XWR template)

Part 3: Model understanding
============================
Amber forcefield was run to investigate the structure: check_ambertools.sh on 3KZ8_DNA_BCD.pdb (modify bash script)

./check_ambertools.sh

Here a range of missing residues and missing heavy atoms was described. 


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
           Running: mod10.1 modeller-multiple3.py

Part 5: Clean up and prepare for analysis
====================================================================================
mkdir models
mv *B9999* models
rm *V9999*
rm *D000* 
