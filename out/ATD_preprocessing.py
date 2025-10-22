#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 23:01:53 2025

@author: val
"""

# If needed:
# %pip install biopython requests

from pathlib import Path
import io, subprocess, shutil, requests
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, PPBuilder, Select
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
from Bio import AlignIO

# --------- paths ---------
DATA = Path("data"); OUT = Path("out")
DATA.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)

TARGET_PDB = Path("out/7EU8_NR1_ATD_only.pdb")   # <-- your target ATD PDB from 7EU8
PDB_8ZH7 = DATA / "8ZH7.pdb"
PDB_5TQ2 = DATA / "5TQ2.pdb"

# Trimmed, single-chain ATDs (templates) MODELLER will read
TEMPL_8ZH7 = OUT / "8ZH7_ATD.pdb"
TEMPL_5TQ2 = OUT / "5TQ2_ATD.pdb"

# Alignment outputs
FASTA_ALL   = OUT / "atd_multi.fasta"
ALN_FASTA   = OUT / "atd_multi_aligned.fasta"
PIR_MULTI   = OUT / "atd_multi.ali"

# --------- utils ---------
def fetch(url, dest: Path):
    if not dest.exists():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        dest.write_bytes(r.content)

def chain_sequence(chain):
    ppb = PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(chain, aa_only=True):
        seq += str(pp.get_sequence())
    return seq

def best_chain_by_identity(pdb_path: Path, target_seq: str, min_len=100):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))
    aligner = PairwiseAligner(); aligner.mode="global"
    try: aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception: pass
    aligner.open_gap_score, aligner.extend_gap_score = -10, -0.5
    best = (None, -1.0, None)  # (chain_id, ident, chain_obj)
    for m in s:
        for ch in m:
            seq = chain_sequence(ch)
            if len(seq) < min_len:
                continue
            aln = aligner.align(seq, target_seq)[0]
            try:
                A,B = str(aln.seqA), str(aln.seqB)
            except AttributeError:
                A,B = map(str, aln[:2])
            nongap = sum(1 for x,y in zip(A,B) if x!='-' and y!='-')
            ident  = sum(1 for x,y in zip(A,B) if x==y and x!='-' and y!='-') / max(nongap,1)
            if ident > best[1]:
                best = (ch.id, ident, ch)
    return best  # (chain_id, identity, chain_obj)

class KeepChainProteinOnly(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    def accept_residue(self, residue):
        hetflag, seqnum, icode = residue.id
        return 1 if hetflag == ' ' else 0
    def accept_atom(self, atom):
        return 1

def write_chain_only(pdb_in: Path, chain_id: str, pdb_out: Path):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_in))
    io = PDBIO(); io.set_structure(s)
    io.save(str(pdb_out), KeepChainProteinOnly(chain_id))

def seq_from_first_chain(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", str(pdb_path))
    for m in s:
        for ch in m:
            seq = chain_sequence(ch)
            if seq:
                return ch.id, seq
    raise RuntimeError(f"No protein sequence found in {pdb_path}")

# --------- target sequence from your 7EU8 ATD PDB ---------
if not TARGET_PDB.exists():
    raise FileNotFoundError(f"Missing target PDB: {TARGET_PDB}")

t_chain, t_seq = seq_from_first_chain(TARGET_PDB)
print(f"[target] chain {t_chain}, length {len(t_seq)}")

# --------- download or use local templates ---------
fetch("https://files.rcsb.org/download/8ZH7.pdb", PDB_8ZH7)
fetch("https://files.rcsb.org/download/5TQ2.pdb", PDB_5TQ2)

# --------- pick best template chain vs target and write chain-only PDBs ---------
ch8, id8, _ = best_chain_by_identity(PDB_8ZH7, t_seq, min_len=100)
ch5, id5, _ = best_chain_by_identity(PDB_5TQ2, t_seq, min_len=100)
print(f"[8ZH7] best chain {ch8}, identity={id8:.1%}")
print(f"[5TQ2] best chain {ch5}, identity={id5:.1%}")

write_chain_only(PDB_8ZH7, ch8, TEMPL_8ZH7)
write_chain_only(PDB_5TQ2, ch5, TEMPL_5TQ2)

# --------- build multi-FASTA: target + two templates ---------
def pdb_chain_seq(pdb_path: Path):
    _, s = seq_from_first_chain(pdb_path)
    return s

seq_8 = pdb_chain_seq(TEMPL_8ZH7)
seq_5 = pdb_chain_seq(TEMPL_5TQ2)

records = [
    SeqIO.SeqRecord(SeqIO.Seq(t_seq), id="target_7EU8_ATD", description=""),
    SeqIO.SeqRecord(SeqIO.Seq(seq_8), id="templ_8ZH7_ATD", description=f"chain {ch8}"),
    SeqIO.SeqRecord(SeqIO.Seq(seq_5), id="templ_5TQ2_ATD", description=f"chain {ch5}"),
]
SeqIO.write(records, FASTA_ALL, "fasta")
print("Wrote", FASTA_ALL)

# --------- align (MAFFT/MUSCLE preferred; fallback to Biopython) ---------
def run_align(in_fa: Path, out_fa: Path):
    if shutil.which("mafft"):
        aln = subprocess.check_output(["mafft", "--auto", str(in_fa)], text=True)
        out_fa.write_text(aln)
        return "mafft"
    if shutil.which("muscle"):
        subprocess.check_call(["muscle", "-align", str(in_fa), "-output", str(out_fa)])
        return "muscle"
    # Fallback: Biopython pairwise progressive (crude)
    recs = list(SeqIO.parse(str(in_fa), "fasta"))
    aligner = PairwiseAligner(); aligner.mode="global"
    try: aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception: pass
    # align templates to target, then merge
    ref = str(recs[0].seq)
    aligned = [(recs[0].id, ref)]
    for r in recs[1:]:
        aln = aligner.align(ref, str(r.seq))[0]
        try: A,B = str(aln.seqA), str(aln.seqB)
        except AttributeError: A,B = map(str, aln[:2])
        aligned.append((r.id, B))
        ref = A  # keep growing reference (simple)
    with out_fa.open("w") as fh:
        for name, s in aligned:
            fh.write(f">{name}\n{s}\n")
    return "biopython_fallback"

tool = run_align(FASTA_ALL, ALN_FASTA)
print("Alignment tool:", tool, "| Wrote", ALN_FASTA)

# --------- write PIR multi-template alignment for MODELLER ---------
aln = AlignIO.read(ALN_FASTA, "fasta")
seqs = {rec.id: str(rec.seq) for rec in aln}

def write_pir(path: Path, seq_target: str, seq_8: str, seq_5: str):
    def fmt(seq):
        return (seq.replace("*","") + "*").replace("\n","")
    with path.open("w") as fh:
        fh.write("C; Multi-template alignment: 8ZH7 + 5TQ2 -> target_7EU8_ATD\n")
        fh.write("C; Gaps are '-' ; '*' ends each sequence for MODELLER\n")
        fh.write("\n")
        fh.write(">P1;8ZH7_ATD\n")
        fh.write("structureX:8ZH7_ATD: : : : : : :\n")
        fh.write(fmt(seq_8) + "\n\n")
        fh.write(">P1;5TQ2_ATD\n")
        fh.write("structureX:5TQ2_ATD: : : : : : :\n")
        fh.write(fmt(seq_5) + "\n\n")
        fh.write(">P1;target_7EU8_ATD\n")
        fh.write("sequence:target_7EU8_ATD: : : : : : :\n")
        fh.write(fmt(seq_target) + "\n")
    print("Wrote PIR:", path)

# Ensure the PDB filenames match the P1; names (8ZH7_ATD.pdb, 5TQ2_ATD.pdb)
write_pir(PIR_MULTI, seqs["target_7EU8_ATD"], seqs["templ_8ZH7_ATD"], seqs["templ_5TQ2_ATD"])
print("Ready for MODELLER: atd_multi.ali + 8ZH7_ATD.pdb + 5TQ2_ATD.pdb")
