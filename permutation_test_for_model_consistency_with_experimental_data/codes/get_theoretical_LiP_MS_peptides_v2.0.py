#!/usr/bin/env python3
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import mdtraj as mdt
import sys,os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle
import time
import itertools

if len(sys.argv) != 8:
    print('[1] path to all-atom reference pdb')
    print('[2] path to join prob pkl file')
    print('[3] outfile_basename')
    print('[4] num iterations')
    print('[5] nproc')
    print('[6] path to PKcutsite observed probability')
    print('[7] path to observed AA prob across proteome')
    quit()

version = 2.0
start_time = time.time()

#get SASA for PDB
traj = mdt.load(sys.argv[1])
sasa = mdt.shrake_rupley(traj, mode='residue').T

with open(sys.argv[2], 'rb') as fh:
    join_prob_dict = pickle.load(fh)
print(join_prob_dict)
outfile_basename = sys.argv[3]

it_reps = int(sys.argv[4])
nproc = int(sys.argv[5])

def gen_theo_peps(it, pairs):

    peptide = np.random.choice(np.arange(len(pairs)), size=1)[0]
    pair = pairs[peptide]
    cutsite = pair[0]
    cutsite_AA = seq[cutsite]
    cutsite_pint = Pint[cutsite_AA]

    sub_seq = "".join(seq[pair[0]+1:pair[1]+1])
    len_sub_seq = len(sub_seq)

    cutsite_str = f'{seq[cutsite]}{cutsite+1}'
    #check for number of internal trypsin cut sites
    num_tryp_sites = 0
    for AA in sub_seq:
        if AA == 'R' or AA == 'K':
            num_tryp_sites+=1

    test_stat = [[len_sub_seq, num_tryp_sites]]
    test_H, test_xedges, test_yedges = np.histogram2d([len_sub_seq], [num_tryp_sites], bins=[join_prob_dict['xedges'], join_prob_dict['yedges']])
    test_P = test_H*join_prob_dict['P']
    test_P = np.sum(test_P)
    di_throw_1 = np.random.uniform(low=0.0, high=1.0, size=None)
    di_throw_2 = np.random.uniform(low=0.0, high=1.0, size=None)
    if di_throw_1 <= cutsite_pint:
        if di_throw_2 <= test_P:
            #print('pair, cutsite, cutsite_AA, cutsite_pint, di_throw_1, sub_seq, len_sub_seq, num_tryp_sites, test_P, di_throw_2')
            #print(pair, cutsite, cutsite_AA, cutsite_pint, di_throw_1, sub_seq, len_sub_seq, num_tryp_sites, test_P, di_throw_2)
            #populate output dictionaries
            return [cutsite_str, cutsite_pint, di_throw_1, sub_seq, len_sub_seq, num_tryp_sites, test_P, di_throw_2]



#get AA sequence and find SASA residue indexs
structure = PDBParser().get_structure('2WW4', sys.argv[1])
ppb=PPBuilder()
global seq
seq = np.asarray(ppb.build_peptides(structure)[0].get_sequence())

#get all possible Trypsin cut sites
Tryp_cutsite_resid = np.where(seq == 'R')[0].astype(int)
Tryp_cutsite_resid = np.hstack((Tryp_cutsite_resid,np.where(seq == 'K')[0].astype(int)))
Tryp_cutsite_AA = seq[Tryp_cutsite_resid]

#get all possible PK cut sites
PK_sasa_resid = np.where(sasa > 0)[0].astype(int)
PK_sasa_resid = [r for r in PK_sasa_resid if r not in Tryp_cutsite_resid]
PK_nonsasa_resid = np.arange(0, len(seq))
PK_nonsasa_resid = [r for r in PK_nonsasa_resid if r not in Tryp_cutsite_resid]

#get propensity for PK to cut at a specific site
#AA_probR1min_C.txt
#AA_probR2hr_C.txt
#AA_probR5min_C.txt
#AA_prob_across_pdb_v1.1.txt
Pobs = sys.argv[6]
Pprot = sys.argv[7]

global Pint
Pint = {}
Pobs = np.loadtxt(Pobs, dtype='O', usecols=(0,2))
Pprot = np.loadtxt(Pprot, dtype='O')
for AA,p in Pobs:
    p_Pprot = Pprot[np.where(Pprot[:,0] == AA)][0][1]

    Pint[AA] = float(p)/float(p_Pprot)

Norm_fact = np.sum(list(Pint.values()))
for AA in Pint.keys():
    Pint[AA] /= Norm_fact


np.savetxt('Pint.txt', [[AA,i] for AA,i in Pint.items()], fmt='%s')
print('SAVED: Pint.txt')

#generate half-tryptic fragments
#Peptide lengths between 6 and 144 amino acid residues were allowed with a peptide mass between 350 and 5000 Da
pairs = itertools.product(PK_nonsasa_resid, Tryp_cutsite_resid)
filtered_pairs = []
for pair in pairs:
    pair = sorted(pair)
    diff = pair[1] - pair[0]
    sub_seq = "".join(seq[pair[0]+1:pair[1]+1])
    len_sub_seq = len(sub_seq)
    if len_sub_seq >= 7 and len_sub_seq <= 50:
        #check molecular weight
        MW = ProteinAnalysis(sub_seq).molecular_weight()
        if MW >= 500 or MW <= 5000:
            filtered_pairs += [pair]
print(len(filtered_pairs))
del pairs

#pk_any_it_outdata_list = list(Parallel(n_jobs=nproc)(delayed(gen_theo_peps)(it, PK_nonsasa_resid, Tryp_cutsite_resid) for it  in range(it_reps)))
pk_any_it_outdata_list = list(Parallel(n_jobs=nproc)(delayed(gen_theo_peps)(it, filtered_pairs) for it  in range(it_reps)))
pk_any_it_outdata_list = [v for v in pk_any_it_outdata_list if v is not None]
out_seq = []
out_cutsites = []
out_len_nRK = []
all_len_nRK = []
out_dict = {}
for data in pk_any_it_outdata_list:
    print(data)
    ss,l,nRK = data[3:6]
    cutsite = data[0]

    all_len_nRK += [[l, nRK]]
    if ss not in out_seq:
        out_seq += [ss]
        out_cutsites += [cutsite]
        out_len_nRK += [[l, nRK]]

out_len_nRK = np.asarray(out_len_nRK)
all_len_nRK = np.asarray(all_len_nRK)

tot_H, xedges, yedges = np.histogram2d(all_len_nRK[:,0], all_len_nRK[:,1], bins=[join_prob_dict['xedges'], join_prob_dict['yedges']])
H, xedges, yedges = np.histogram2d(out_len_nRK[:,0], out_len_nRK[:,1], bins=[join_prob_dict['xedges'], join_prob_dict['yedges']])
print('total_theoretical join dist')
print(tot_H/np.sum(tot_H))
print('theoretical join dist')
print(H/np.sum(H))
print('obs join dist')
print(join_prob_dict['P'])
out_dict['any_pk_site'] = {}
out_dict['any_pk_site']['cutsites'] = out_cutsites
out_dict['any_pk_site']['sub_seqs'] = out_seq
out_dict['any_pk_site']['len_vs_nRK'] = out_len_nRK
out_dict['any_pk_site']['len_vs_nRK_Hist'] = H
out_dict['any_pk_site']['len_vs_nRK_Prob'] = H/np.sum(H)
out_dict['any_pk_site']['tot_len_vs_nRK_Hist'] = tot_H
out_dict['any_pk_site']['tot_len_vs_nRK_Prob'] = tot_H/np.sum(tot_H)

with open(f'{outfile_basename}_any_pk_site_data.txt', 'w') as fh:
    fh.write('#cutsite\n')
    for rd in out_cutsites:
        fh.write(f'{rd}\n')
print(f'SAVED: {outfile_basename}_any_pk_site_data.txt')


with open(f'{outfile_basename}_theoretical_LiP_MS_peptides_output.pkl', 'wb') as fh:
    pickle.dump(out_dict, fh)
print(f'SAVED: {outfile_basename}_theoretical_LiP_MS_peptides_output.pkl')
print(f'NORMAL TERMINATION @ {time.time() - start_time}')
