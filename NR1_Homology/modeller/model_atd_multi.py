
from modeller import * 
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller import soap_protein_od

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', '/Users/val/Projects/NMDAR/NR1_Homology/modeller/out', '/Users/val/Projects/NMDAR/NR1_Homology/modeller/out/../data']

# --- define class ---
class MyModel(automodel):
    def select_atoms(self):
        return selection(self.chains[0])

a = MyModel(env,
            alnfile=r'/Users/val/Projects/NMDAR/NR1_Homology/modeller/out/atd_multi.ali',
            knowns=('8ZH7_ATD','5TQ2_ATD'),
            sequence='target_7EU8_ATD')

# Send all outputs explicitly to your modeller/out directory
a.outputs_dir = r'/Users/val/Projects/NMDAR/NR1_Homology/modeller/out'
a.starting_model = 1
a.ending_model   = 5
a.md_level = refine.very_fast

# Ask MODELLER to compute DOPE and GA341
a.assess_methods = (assess.DOPE, assess.GA341)

a.make()

# --- summary block (robust to both DOPE / MolPDF) ---
ok = [m for m in a.outputs if m['failure'] is None]

def get_score(m):
    if 'DOPE score' in m:
        return ('DOPE', m['DOPE score'])
    if 'GA341 score' in m:
        return ('GA341', m['GA341 score'])
    if 'MolPDF' in m:
        return ('MolPDF', m['MolPDF'])
    if 'molpdf' in m:
        return ('molpdf', m['molpdf'])
    return ('unknown', 1e99)

ok.sort(key=lambda m: get_score(m)[1])

for m in ok:
    label, val = get_score(m)
    print("MODEL %s  %s=%.3f" % (m['name'], label, val))
