
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', './out', './data']

# Templates must match the P1; names in out/atd_multi.ali
class MyModel(automodel):
    def select_atoms(self):
        # refine whole model (fast); you can target loop regions if needed
        return selection(self.chains[0])

a = MyModel(env,
            alnfile='out/atd_multi.ali',
            knowns=('8ZH7_ATD','5TQ2_ATD'),
            sequence='target_7EU8_ATD')

a.starting_model = 1
a.ending_model   = 5         # build 5 models
a.md_level = refine.very_fast
a.make()

# Print a DOPE-sorted summary
ok = [m for m in a.outputs if m['failure'] is None]
ok.sort(key=lambda x: x['DOPE score'])
for m in ok:
    print(f"MODEL {m['name']}  DOPE={m['DOPE score']:.1f}")
