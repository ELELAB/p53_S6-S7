#comparative modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

class MyModel(automodel):
      def user_after_single_model(self):
        # each model:
        self.rename_segments(segment_ids=('B','B','C','D'), renumber_residues=[91, 300, 1, 1])
      def select_atoms(self):
          return selection(self.residue_range('1:A', '6:A'), self.residue_range('95:A', '97:A'),self.atoms['NZ:11:A'],self.atoms['CE:11:A'],self.atoms['CG:14:A'],self.atoms['CD:14:A'],self.atoms['NE2:14:A'],self.atoms['OE1:14:A'],self.atoms['CD:20:A'],self.atoms['NE:20:A'],self.atoms['CZ:20:A'],self.atoms['NH1:20:A'],self.atoms['NH2:20:A'],self.atoms['CG:58:A'],self.atoms['OD1:58:A'],self.atoms['OD2:58:A'],self.atoms['CG:98:A'],self.atoms['CD1:98:A'],self.atoms['CD2:98:A'],self.atoms['CG:111:A'],self.atoms['CD1:111:A'],self.atoms['CD2:111:A'],self.atoms['NE:112:A'],self.atoms['CZ:112:A'],self.atoms['NH1:112:A'],self.atoms['NH2:112:A'],self.atoms['CD:114:A'],self.atoms['OE1:114:A'],self.atoms['OE2:114:A'],self.atoms['CG:131:A'],self.atoms['CD:131:A'],self.atoms['OE1:131:A'],self.atoms['OE2:131:A'],self.atoms['CG:134:A'],self.atoms['CD:134:A'],self.atoms['OE1:134:A'],self.atoms['OE2:134:A'],self.atoms['CG1:135:A'],self.atoms['CG2:135:A'],self.atoms['CG:138:A'],self.atoms['OD1:138:A'],self.atoms['OD2:138:A'],self.atoms['CD:197:A'],self.atoms['OE1:197:A'],self.atoms['OE2:197:A'])
      def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Restrain the specified distance to X angstroms (st. dev.=0.1)
#       Use a harmonic potential and X-Y distance group.
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:84:A'],
                                                         at['CG:1:A']),
                               mean=4.0, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:84:A'],
                                                         at['CH2:1:A']),
                               mean=4.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:121:A'],
                                                         at['OG:4:A']),
                               mean=2.7, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['OG:4:A'],
                                                         at['N:6:A']),
                               mean=3.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:80:A'],
                                                         at['N:4:A']),
                               mean=2.8, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CB:2:A'],
                                                         at['CG:122:A']),
                               mean=5.0, stdev=0.1))
env = environ(rand_seed=-1289)
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

env.io.hetatm = True

a = MyModel(env,
            alnfile  = 'xxx.ali', # alignment filename
            knowns   = ('3KZ8_DNA_BCD', '2XWRa_aligned', '2XWRa_aligned_2'),     # codes of the templates
            sequence = 'p53_3KZ8_DNA_BCD_2xwra_91-289')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 100               # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do the actual comparative modeling
