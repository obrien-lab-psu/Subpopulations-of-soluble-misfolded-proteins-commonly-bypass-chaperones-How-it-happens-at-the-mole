#!/usr/bin/env python3
import pandas as pd
import numpy as np
import mdtraj as mdt
import sys,os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle

if len(sys.argv) != 5:
    print('[1] path to all-atom reference pdb')
    print('[2] outfile must end in .pkl')
    print('[3] path to dcd')
    print('[4] path to dcd topology')
    quit()

#get SASA for PDB
ref_traj = mdt.load(sys.argv[1])
ref_sasa = mdt.shrake_rupley(ref_traj, mode='residue').T

msm_traj = mdt.load(sys.argv[3], top=sys.argv[4])
#msm_traj = mdt.load(sys.argv[3])
msm_sasa = mdt.shrake_rupley(msm_traj, mode='residue').T
print(ref_sasa.shape)
print(ref_sasa[:10])
print(msm_sasa.shape)
print(msm_sasa[:10])
diff = msm_sasa - ref_sasa
perc_diff = 100*((msm_sasa - ref_sasa)/ref_sasa)
print(diff.shape)
print(diff[:10])
diff = diff.T
msm_sasa = msm_sasa.T
ref_sasa = ref_sasa.T

#get AA sequence and find SASA residue indexs
structure = PDBParser().get_structure('2WW4', sys.argv[1])
ppb=PPBuilder()
seq = np.asarray(ppb.build_peptides(structure)[0].get_sequence())
print(seq)

outdict = {'seq':seq, 'ref_sasa':ref_sasa, 'msm_sasa':msm_sasa, 'diff':diff, 'perc_diff':perc_diff}
with open(sys.argv[2], 'wb') as fh:
    pickle.dump(outdict, fh)

print(f'SAVED: {sys.argv[2]}')

