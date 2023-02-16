#!/usr/bin/env python3
import sys,os
import numpy as np
import glob

if len(sys.argv) != 2:
    print('[1] path to file containing trapped traj info')
    quit()

#trapped_info = {x[0]:np.asarray(x[1].split(',')).astype(int) for x in np.loadtxt(sys.argv[1], dtype='O', delimiter=',')}
trapped_info = {}
for line in open(sys.argv[1]).readlines()[1:]:
    line = line.split(',')
    pdb = line[0]
    traj = int(line[1])

    if pdb in ['3hwo', '2fym', '4a2c', '1k7j', '1p7l', '1a69']:
        if pdb not in trapped_info:
            trapped_info[pdb] = [traj]
        else:
            trapped_info[pdb] += [traj]

for k,v in trapped_info.items():
    print(k,v)

for pdb, trajs in trapped_info.items():
    trapped = 0
    trapped_and_change_ent = 0
    trajs  = np.asarray(trajs)
    for t in trajs:
        trapped += 1
        t_Qfile = f'entanglement_analysis_2.0/output/p{pdb}_post_trans_t{t}_GQ.txt'

        G = np.loadtxt(t_Qfile, usecols=(2))

        num_nonzero = np.nonzero(G)[0].size

        if num_nonzero/G.size >= 0.5:
            trapped_and_change_ent += 1

    print(f'{pdb} | %misfolded: {trapped/50} | %changeEnt: {trapped_and_change_ent/trajs.size}')
