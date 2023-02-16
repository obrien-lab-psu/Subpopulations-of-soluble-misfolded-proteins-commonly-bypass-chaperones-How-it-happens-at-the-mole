#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys,os

if len(sys.argv) != 3:
    print('[1] = input crd')
    print('[2] = output crd')
    quit()

u = mda.Universe(sys.argv[1])
print(u)
protein = u.select_atoms('all')

for i,res in enumerate(protein.residues):
    print(i+1, res, res.resid)
    res.resid = i + 1
    print(i+1, res, res.resid)
    res.segment.segid= 'A'

with mda.Writer(sys.argv[2], protein.n_atoms) as W:
    W.write(protein)

quit()
