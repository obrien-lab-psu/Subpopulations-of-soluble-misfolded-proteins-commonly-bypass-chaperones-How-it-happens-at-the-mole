#!/usr/bin/env python3

import sys,os
import numpy as np
import MDAnalysis as mda

if len(sys.argv) != 6:
    print('[1] = dcd')
    print('[2] = psf')
    print('[3] = working dir')
    print('[4] = outfile basename')
    print('[5] = AA PDB')
    sys.exit()

dcd = sys.argv[1]
psf = sys.argv[2]
wd = sys.argv[3]
bn = sys.argv[4]
aa_pdb = sys.argv[5]
#print(dcd, psf, wd, bn, aa_pdb)

if not os.path.exists(wd):
    os.mkdir(wd)
    #print(f'Made: {wd}')

u = mda.Universe(psf,dcd)
#print(u)

protein = u.select_atoms('all')
for ts in u.trajectory:

    with mda.Writer(f'{wd}{bn}_{ts.frame}.pdb', protein.n_atoms) as W:
        W.write(protein)
    #print(f'Saved: {wd}{bn}_{ts.frame}.pdb')

    backmap_output = f'python codes/backmap_pulchra_only.py -i {aa_pdb} -c {wd}{bn}_{ts.frame}.pdb'
    print(backmap_output)


