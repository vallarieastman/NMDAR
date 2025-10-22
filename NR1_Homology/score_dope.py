#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 13:55:50 2025

@author: val
"""

from collections import Counter, defaultdict

#DSSP  = "/Users/val/Downloads/scrambled_NR1_ESM/scrambled_NR1_ESMfold.dssp"
pdb = "/Users/val/Downloads/Vallari_Eastman/model_02.pdb"

from modeller import environ, model
from modeller.scripts import complete_pdb
from modeller import assessment

def dope_score(pdbfile):
    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    mdl = complete_pdb(env, pdbfile)
    atmsel = selection(mdl)
    s = atmsel.assess_dope()
    return s

for pdb in [pdb]:
    print(pdb, "DOPE:", dope_score(pdb))



