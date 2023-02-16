#!/bin/usr/python
################################################################################################################
#TITLE: Test MDAnalysis RMSD FIT of each frame in ref traj to every frame in another traj
#
#AUTHOR: IAN S 2019.12.17
################################################################################################################

import sys
import numpy as np
import re
import time
from scipy import stats
import MDAnalysis
from scipy.spatial.distance import cdist

#check for arguments
if  len(sys.argv)!=5:
    print('[1] = psf')
    print('[2] = dcd')
    print('[3] = outfile')
    print('[4] = client chain ID')
    quit()

PSF = sys.argv[1]
DCD = sys.argv[2]
outfile = sys.argv[3]
client_chains = sys.argv[4]
print(PSF, DCD, outfile)

#load universe with dcd and psf
u = MDAnalysis.Universe(PSF,DCD)

#correct the chain ID element of the universe
chains = []
for res in u.residues:
    chain = res.resname[0]
    chains.append(chain)
    res.segment.segid = chain

#sepearate client and nonclient
client = u.select_atoms(f'segid {client_chains}')
nonclient = u.select_atoms(f'not segid {client_chains}')
print(f'CLIENT CHAINS: {client}')
print(f'NONCLIENT CHAINS: {nonclient}')

#loop through traj and calculate the binding events
outdata = []
print(f'TOTAL FRAMES TO ANALYZE: {len(u.trajectory)}')
for ts in u.trajectory:

    #get coordinates
    client_coor = client.positions
    nonclient_coor = nonclient.positions
    #print(client_coor.shape)
    #print(nonclient_coor.shape)

    #get coordinate distances
    dist = cdist(client_coor, nonclient_coor)
    #print(dist.shape)

    #determine if there are any distances less than 8A
    contacts = np.where(dist <= 8)
    num_contacts = len(contacts[0])

    if num_contacts >= 150:
        num_contacts = 1
    else:
        num_contacts = 0

    outdata.append([ts.frame, num_contacts])

outdata = np.vstack(outdata)
np.savetxt(outfile,outdata)
print(f'SAVED: {outfile}')
###################################################################################################################
