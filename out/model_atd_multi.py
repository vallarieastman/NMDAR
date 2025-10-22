#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 22:57:58 2025

@author: val
"""

from pathlib import Path

script = r"""
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', './out', './data']

# Make sure the PDB basenames match the P1; names in the PIR file:
#   8ZH7_ATD.pdb, 5TQ2_ATD.pdb  (in ./out)
# And the alignment file is: atd_multi.ali

class MyLoopRefine(automodel):
    def select_atoms(self):
        # You can refine loops by returning a selection here; default = whole model
        return selection(self.chains[0])

# Known templates must match P1; entry IDs in the PIR
a = MyLoopRefine(env,
                 alnfile='out/atd_multi.ali',
                 knowns=('8ZH7_ATD','5TQ2_ATD'),
                 sequence='target_7EU8_ATD')

a.starting_model = 1
a.ending_model   = 5        # build 5 models
a.md_level = refine.very_fast  # quick loop refinement
a.make()

# Write a summary table sorted by DOPE
# Write a summary table sorted by DOPE
ok = [m for m in a.outputs if m['failure'] is None]
ok.sort(key=lambda x: x['DOPE score'])
for m in ok:
    print("MODEL %s  DOPE=%.1f" % (m['name'], m['DOPE score']))
"""

Path("model_atd_multi.py").write_text(script)
print("Wrote model_atd_multi.py")
print("Run with:  mod10.7 model_atd_multi.py")
