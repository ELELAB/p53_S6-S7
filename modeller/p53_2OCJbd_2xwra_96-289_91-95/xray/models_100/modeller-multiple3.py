#comparative modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

class MyModel(automodel):
      def user_after_single_model(self):
        # each model:
        self.rename_segments(segment_ids=('B', 'D'), renumber_residues=[91, 91])
      def select_atoms(self):
        return selection(self.residue_range('1:A', '6:A'), self.residue_range('201:B', '206:B'))
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
#second chain
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:284:B'],
                                                         at['CG:201:B']),
                               mean=4.0, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:284:B'],
                                                         at['CH2:201:B']),
                               mean=4.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:321:B'],
                                                         at['OG:204:B']),
                               mean=2.7, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['OG:204:B'],
                                                         at['N:206:B']),
                               mean=3.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:280:B'],
                                                         at['N:204:B']),
                               mean=2.8, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CB:202:B'],
                                                         at['CG:322:B']),
                               mean=5.0, stdev=0.1))

env = environ(rand_seed=-1208)
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

env.io.hetatm = True

a = MyModel(env,
            alnfile  = 'xxx.ali', # alignment filename
            knowns   = ('2OCJ_clean_BD', '2XWR_clean_A_alignB', '2XWR_clean_A_alignD'),     # codes of the templates
            sequence = 'p53_2OCJbd_2xwra_91-289')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 100               # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do the actual comparative modeling
