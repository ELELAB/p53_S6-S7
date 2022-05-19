#comparative modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

class MyModel(automodel):
      def user_after_single_model(self):
        # each model:
        self.rename_segments(segment_ids=('A','A','B','B'), renumber_residues=[91, 300, 91, 300])
      def select_atoms(self):
          return selection(self.residue_range('1:A', '6:A'), self.residue_range('95:A', '97:A'),
                  self.atoms['NZ:11:A'],self.atoms['CE:11:A'],
                  self.atoms['CG:14:A'],self.atoms['CD:14:A'],self.atoms['NE2:14:A'],self.atoms['OE1:14:A'],
                  self.atoms['CD:20:A'],self.atoms['NE:20:A'],self.atoms['CZ:20:A'],self.atoms['NH1:20:A'],self.atoms['NH2:20:A'],
                  self.atoms['CG:58:A'],self.atoms['OD1:58:A'],self.atoms['OD2:58:A'],
                  self.atoms['CG:98:A'],self.atoms['CD1:98:A'],self.atoms['CD2:98:A'],
                  self.atoms['CG:111:A'],self.atoms['CD1:111:A'],self.atoms['CD2:111:A'],
                  self.atoms['NE:112:A'],self.atoms['CZ:112:A'],self.atoms['NH1:112:A'],self.atoms['NH2:112:A'],
                  self.atoms['CD:114:A'],self.atoms['OE1:114:A'],self.atoms['OE2:114:A'],
                  self.atoms['CG:131:A'],self.atoms['CD:131:A'],self.atoms['OE1:131:A'],self.atoms['OE2:131:A'],
                  self.atoms['CG:134:A'],self.atoms['CD:134:A'],self.atoms['OE1:134:A'],self.atoms['OE2:134:A'],
                  self.atoms['CG1:135:A'],self.atoms['CG2:135:A'],
                  self.atoms['CG:138:A'],self.atoms['OD1:138:A'],self.atoms['OD2:138:A'],
                  self.atoms['CD:197:A'],self.atoms['OE1:197:A'],self.atoms['OE2:197:A'],
                  self.residue_range('201:C', '206:C'), 
                  self.residue_range('295:C', '297:C'),
                  self.atoms['NZ:211:C'],self.atoms['CE:211:C'],
                  self.atoms['CG:214:C'],self.atoms['CD:214:C'],self.atoms['NE2:214:C'],self.atoms['OE1:214:C'],
                  self.atoms['CD:220:C'],self.atoms['NE:220:C'],self.atoms['CZ:220:C'],self.atoms['NH1:220:C'],self.atoms['NH2:220:C'],
                  self.atoms['CG:258:C'],self.atoms['OD1:258:C'],self.atoms['OD2:258:C'],
                  self.atoms['CG:298:C'],self.atoms['CD1:298:C'],self.atoms['CD2:298:C'],
                  self.atoms['CG:311:C'],self.atoms['CD1:311:C'],self.atoms['CD2:311:C'],
                  self.atoms['NE:312:C'],self.atoms['CZ:312:C'],self.atoms['NH1:312:C'],self.atoms['NH2:312:C'],
                  self.atoms['CD:314:C'],self.atoms['OE1:314:C'],self.atoms['OE2:314:C'],
                  self.atoms['CG:331:C'],self.atoms['CD:331:C'],self.atoms['OE1:331:C'],self.atoms['OE2:331:C'],
                  self.atoms['CG:334:C'],self.atoms['CD:334:C'],self.atoms['OE1:334:C'],self.atoms['OE2:334:C'],
                  self.atoms['CG1:335:C'],self.atoms['CG2:335:C'],
                  self.atoms['CG:338:C'],self.atoms['OD1:338:C'],self.atoms['OD2:338:C'],
                  self.atoms['CD:397:C'],self.atoms['OE1:397:C'],self.atoms['OE2:397:C'])
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
#       Restrain the specified distance to X angstroms (st. dev.=0.1)
#       Use a harmonic potential and X-Y distance group.
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:284:C'],
                                                         at['CG:201:C']),
                               mean=4.0, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CZ:284:C'],
                                                         at['CH2:201:C']),
                               mean=4.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:321:C'],
                                                         at['OG:204:C']),
                               mean=2.7, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['OG:204:C'],
                                                         at['N:206:C']),
                               mean=3.3, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['O:280:C'],
                                                         at['N:204:C']),
                               mean=2.8, stdev=0.1))
        rsr.add(forms.Gaussian(group=physical.xy_distance,
                               feature=features.Distance(at['CB:202:C'],
                                                         at['CG:322:C']),
                               mean=5.0, stdev=0.1))
        
env = environ(rand_seed=-1056)
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

env.io.hetatm = True

a = MyModel(env,
            alnfile  = 'xxx.ali', # alignment filename
            knowns   = ('3kz8_01', '2XWRa_aligned_A', '2XWRa_aligned_B'),     # codes of the templates
            sequence = 'p53_3KZ8_AB_2xwra_91-289')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 100               # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do the actual comparative modeling
