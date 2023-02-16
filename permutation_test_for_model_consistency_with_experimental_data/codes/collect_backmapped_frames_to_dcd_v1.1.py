#!/usr/bin/env python3
import MDAnalysis as mda
import sys,os
import glob

if len(sys.argv) != 3:
    print('[1] path to inputt backmapped pdbs files')
    print('[2] path to outfile')
    quit()

pdbs = glob.glob(sys.argv[1])
print(pdbs)
pdbs = sorted(pdbs, key=lambda x:int(x.split('/')[-1].split('_')[1]))
print(pdbs)
outfile = sys.argv[2]
print(outfile)

u = mda.Universe(pdbs[0], pdbs)
print(u)

protein = u.select_atoms('all')
with mda.Writer(f'{outfile}', protein.n_atoms) as W:
    for ts in u.trajectory:

        W.write(protein)
