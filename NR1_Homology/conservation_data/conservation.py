#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 20:49:27 2025

@author: val
"""

from pymol import cmd
import csv

def map_conservation(csv_path,
                     obj="NR1",
                     rescol="human_resnum",
                     scorecol="conservation_score_0to9",
                     offset=0):
    """
    csv_path : path to CSV with per-residue conservation
    obj      : object name in PyMOL
    rescol   : column with human residue number (int-like)
    scorecol : column with conservation score (we'll map this -> b-factor)
    offset   : add this to human residue number to match PDB numbering
    """

    # make sure offset is an int, even if PyMOL passed a string
    try:
        offset = int(offset)
    except Exception:
        print(f"[map_conservation] Warning: couldn't cast offset={offset} to int, using 0")
        offset = 0

    # 1. collect which residue numbers actually exist in the object
    pdb_resis = []
    cmd.iterate(f"{obj} and name CA", "pdb_resis.append(int(resi))", space={"pdb_resis": pdb_resis})
    pdb_resis_set = set(pdb_resis)
    if len(pdb_resis_set) == 0:
        print(f"[map_conservation] ERROR: object '{obj}' has no residues (check obj name).")
        return

    # helpful for debugging:
    print(f"[map_conservation] PDB residue index range for {obj}: {min(pdb_resis_set)}–{max(pdb_resis_set)} (CA atoms)")

    # 2. read CSV and assign b-factors only when residue exists in PDB
    missing_count = 0
    assigned_count = 0

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row.get(rescol) or not row.get(scorecol):
                continue

            resi_human = int(float(row[rescol]))
            score      = float(row[scorecol])

            resi_pdb   = resi_human + offset

            if resi_pdb in pdb_resis_set:
                sel = f"{obj} and resi {resi_pdb}"
                cmd.alter(sel, f"b={score}")
                assigned_count += 1
            else:
                missing_count += 1

    cmd.rebuild()

    print(f"[map_conservation] Assigned b-factors for {assigned_count} residues.")
    if missing_count > 0:
        print(f"[map_conservation] Skipped {missing_count} residues not present in {obj} (e.g. outside ATD).")

cmd.extend("map_conservation", map_conservation)

