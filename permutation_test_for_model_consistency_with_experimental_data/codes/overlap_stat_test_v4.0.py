#!/usr/bin/env python3
#calculate the pvalue for the probability that an entangled state i has greater overlap at a given time point with the set of sig peptides than in the observed state.
import numpy as np
from joblib import Parallel, delayed
import sys,os
import glob
import pickle
import pandas as pd
import re
from scipy.spatial import distance
from itertools import product
import time
import ast
from statsmodels.stats.multitest import fdrcorrection

np.set_printoptions(precision=3)
pd.options.display.max_columns = 2000
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 2000)
pd.set_option('display.width', 500)

if len(sys.argv) != 10:
    print('[1] path to ent  state .pkl file')
    print('[2] path to LiPMS_peptide file')
    print('[3] path to outfile')
    print('[4] buffer')
    print('[5] path to theoretical peptide file')
    print('[6] path to sasa change file')
    print('[7] max_res')
    print('[8] reps')
    print('[9] number processors')
    sys.exit()

def load_rep_ent_data(rep_ent):
    print('--------------------------------------------------------------------------------------------------------------')
    print('\nLoading all rep ent')
    rep_ent_data = {}
    for row_i, row in enumerate(rep_ent):
        for nc, ent_info in row.items():

            surr = np.hstack(ent_info[1][3])
            frame = ent_info[0]

            ent_res = set(surr)
            nc_res = set(np.hstack((np.arange(nc[0]-res_buffer, nc[0]+res_buffer), np.arange(nc[1]-res_buffer, nc[1]+res_buffer))))
            ent_res = ent_res.union(nc_res)

            rep_ent_data[row_i] = {}
            rep_ent_data[row_i]['info'] = nc
            rep_ent_data[row_i]['res'] = ent_res
            rep_ent_data[row_i]['frame'] = frame

    return rep_ent_data

def jaccard(A,B):
    return len(A.intersection(B))/len(A.union(B))

def load_LipMS_data(file_path):
    df = pd.read_csv(file_path)
    print(df)
    LipMS_change_data = {}
    LipMS_sig_data = {}
    u_timepoints = np.unique(df['Timepoint'])
    sig_peptides = df[df['Significant'] ==True]['proteinaseKsite']
    sig = df[df['Significant'] ==True]

    cutsites = [int(re.sub('\D','',x)) for x in sig['proteinaseKsite']]
    sig['cutsite_idx'] = cutsites
    print(sig.sort_values(by=['cutsite_idx']))

    sig_pk_cutsites, num_sig_peptides = np.unique(sig_peptides, return_counts=True)
    timepoint_2_num_sig_peptides=  {}

    ##########################################################
    #0 Unnamed: 0                                      9636
    #Accession                                     P0ABP8
    #Peptide Sequence    [M].ATPHINAEMGDFADVVLMPGDPLR.[A]
    #proteinaseKsite                                   A2
    #PeptideRatio1                              -0.353492
    #PeptidePValue1                              0.292641
    #Significant                                    False
    #Buffer                                             B
    #Timepoint                                          1
    ##########################################################

    for index, row in df.iterrows():
        peptide = row['proteinaseKsite']
        timepoint = row['Timepoint']
        if timepoint not in timepoint_2_num_sig_peptides:
            timepoint_2_num_sig_peptides[timepoint] = 0

        if '[' not in peptide:
            peptide_range = int(re.sub('\D', '', peptide))
            peptide_range = set(list(np.arange(peptide_range-res_buffer, peptide_range+res_buffer+1)))
            peptide_range = [x for x in peptide_range if x >= 1]
            peptide_range = [x for x in peptide_range if x < max_res]
            #LipMS_all_data[peptide] = peptide_range

            if row['Significant'] == True:
                timepoint_2_num_sig_peptides[timepoint] += 1

                if peptide not in LipMS_sig_data:
                    LipMS_sig_data[peptide] = {}
                    LipMS_sig_data[peptide]['peptide_range'] = []
                    LipMS_sig_data[peptide]['timepoints'] = []
                    LipMS_sig_data[peptide]['qual_change'] = []

                #LipMS_sig_data[timepoint][peptide] = peptide_range
                LipMS_sig_data[peptide]['peptide_range'] += [peptide_range]
                LipMS_sig_data[peptide]['timepoints'] += [timepoint]

                if row['PeptideRatio1'] < 0:
                    LipMS_sig_data[peptide]['qual_change'] += [-1]
                elif row['PeptideRatio1'] > 0:
                    LipMS_sig_data[peptide]['qual_change'] += [1]
                elif row['PeptideRatio1'] == 0:
                    LipMS_sig_data[peptide]['qual_change'] += [0]

    return LipMS_sig_data, u_timepoints, timepoint_2_num_sig_peptides

def load_theo_peps(file_path):
    LipMS_all_data = {}
    data = np.loadtxt(file_path, dtype='O')
    for peptide in data:
        peptide_range = int(re.sub('\D', '', peptide))
        peptide_range = set(list(np.arange(peptide_range-res_buffer, peptide_range+res_buffer+1)))
        peptide_range = [x for x in peptide_range if x >= 1]
        peptide_range = [x for x in peptide_range if x < max_res]
        LipMS_all_data[peptide] = peptide_range

    return LipMS_all_data

def consistency(lip_data, u_timepoints, rep_ent_data, mode):

    if mode == 'T':
        rand_peps = np.random.choice(list(lip_data.keys()), size=num_u_LipMS_sig_peps, replace=False)
        #print(rand_peps)
        num_u_rep_peps = 0
        while num_u_rep_peps != num_u_LipMS_sig_peps:
            rep_LipMS_data = {}
            rep_timepoint_peps = {}
            for timepoint in u_timepoints:
                rep_timepoint_peps[timepoint] = np.random.choice(rand_peps, size=timepoint_2_num_sig_peptides[timepoint], replace=False)

                for peptide in rep_timepoint_peps[timepoint]:
                    if peptide not in rep_LipMS_data:
                        rep_LipMS_data[peptide] = {'peptide_range':[], 'timepoints':[], 'qual_change':[]}

                    rep_LipMS_data[peptide]['peptide_range'] += [lip_data[peptide]]
                    #print(peptide, LipMS_all_data[peptide])
                    rep_LipMS_data[peptide]['timepoints'] += [timepoint]
                    rep_LipMS_data[peptide]['qual_change'] += [np.random.choice([1,-1], size=1)[0]]


            num_u_rep_peps = len(np.unique(np.hstack(list(rep_timepoint_peps.values()))))

        LipMS_sig_data = rep_LipMS_data
        del rep_LipMS_data
        del lip_data

    elif mode == 'O':
        LipMS_sig_data = lip_data
        del lip_data

    total_overlap_consistency = []
    total_exposure_consistency = []
    for obs_sig_pk_cutsite_idx, obs_sig_pk_cutsite in enumerate(LipMS_sig_data.keys()):
        #print('\n---------------------------------------------------------------------------')
        #print(obs_sig_pk_cutsite_idx, obs_sig_pk_cutsite)
        #0 A227
        #peptide_range [[224, 225, 226, 227, 228, 229, 230, 231, 232, 222, 223], [224, 225, 226, 227, 228, 229, 230, 231, 232, 222, 223]]
        #timepoints [1, 5]
        #qual_change [-1, -1]

        #########
        ### get exposure consistency

        #########
        ###get overlap of peptide with each entanglement
        overlap_consistency = []
        exposure_consistency = []
        lip_time_data = LipMS_sig_data[obs_sig_pk_cutsite]['timepoints']
        lip_change_data = LipMS_sig_data[obs_sig_pk_cutsite]['qual_change']
        for timepoint_idx, timepoint in enumerate(u_timepoints):
            #print(timepoint_idx, timepoint)
            if timepoint not in lip_time_data:
                ent_overlap = np.zeros(len(rep_ent_data))
                ent_exposure = np.zeros(len(rep_ent_data))
                overlap_consistency += [ent_overlap[:,None]]
                exposure_consistency += [ent_exposure[:,None]]
                continue

            else:
                ent_overlap = np.zeros(len(rep_ent_data))
                ent_exposure = np.zeros(len(rep_ent_data))

                local_timepoint_idx = np.where(lip_time_data==timepoint)[0][0]
                local_change_data = lip_change_data[local_timepoint_idx]
                lip_res = LipMS_sig_data[obs_sig_pk_cutsite]['peptide_range'][local_timepoint_idx]
                lip_res = np.asarray(lip_res).astype(int)
                for ent_id in rep_ent_data.keys():
                    #print('\n',ent_id)
                    #print(rep_ent_data[ent_id])
                    #print(lip_res)
                    cutsite = np.median(lip_res).astype(int)
                    #print(cutsite)

                    ent_overlap[ent_id] = jaccard(set(rep_ent_data[ent_id]['res']), set(lip_res))

                    #ent_pk_sasa_change = np.mean(sasa_diff[rep_ent_data[ent_id]['frame']][lip_res])
                    ent_pk_sasa_change = np.mean(sasa_diff[rep_ent_data[ent_id]['frame']][cutsite])

                    if ent_pk_sasa_change > 0 and local_change_data > 0:
                        ent_exposure[ent_id] = 1
                    if ent_pk_sasa_change < 0 and local_change_data < 0:
                        ent_exposure[ent_id] = 1

                overlap_consistency += [ent_overlap[:,None]]
                exposure_consistency += [ent_exposure[:,None]]

        overlap_consistency = np.hstack(overlap_consistency)
        total_overlap_consistency += [overlap_consistency]
        exposure_consistency = np.hstack(exposure_consistency)
        total_exposure_consistency += [exposure_consistency]

    #total_overlap_consistency shape = (#LiPMS_PK, #change_ent, #timepoints)
    total_overlap_consistency = np.stack(total_overlap_consistency)
    total_exposure_consistency = np.stack(total_exposure_consistency)
    total_overlap_consistency = np.where(total_overlap_consistency>0, 1, 0)
    #print(total_overlap_consistency.shape, total_exposure_consistency.shape)
    total_obs_overlap_mean_consistency = np.sum(total_overlap_consistency,  axis=(0,1))
    total_obs_exposure_mean_consistency = np.sum(total_exposure_consistency,  axis=(0,1))

    if mode == 'O':
        return total_obs_overlap_mean_consistency, total_obs_exposure_mean_consistency

    elif mode == 'T':

        #print(total_obs_overlap_mean_consistency, total_obs_exposure_mean_consistency)
        #print(obs_overlap_consistency, obs_exposure_consistency)
        condition1 = total_obs_overlap_mean_consistency > obs_overlap_consistency
        condition2 = total_obs_exposure_mean_consistency > obs_exposure_consistency

        found_overlap = 0
        found_overlap_at_max = False
        for timepoint_idx, timepoint in enumerate(u_timepoints):
            if condition1[timepoint_idx] == True:
                found_overlap += 1
                if timepoint == max_timepoint:
                    found_overlap_at_max = True

        #print(found_overlap, found_overlap_at_max)
        if found_overlap_at_max == True and found_overlap >= 2:
            return 1
        else:
            return 0


############
### MAIN ###
############
reps = int(sys.argv[8])
global max_res
max_res = int(sys.argv[7])
start_time = time.time()
outfile = sys.argv[3]
global res_buffer
res_buffer = int(sys.argv[4])
nproc = int(sys.argv[9])

#load SASA change data
# if misfolded is more exposed +1 if native is more exposed -1
#sasa_diff = treatment data - native data so if sasa_diff > 0 then misfolded was more exposed than native
with open(sys.argv[6], 'rb') as fh:
    sasa_change_data = pickle.load(fh)

global sasa_diff
sasa_diff = sasa_change_data['diff']
print(sasa_diff.shape)

#load lipms data
LipMS_sig_data, u_timepoints, timepoint_2_num_sig_peptides = load_LipMS_data(sys.argv[2])
for k,v in LipMS_sig_data.items():
    print(k,v)

num_u_LipMS_sig_peps = len(LipMS_sig_data.keys())
global max_timepoint
max_timepoint = max(u_timepoints)

#load null peptide set
LipMS_all_data = load_theo_peps(sys.argv[5])
for peptide, lipdata in LipMS_sig_data.items():
    print(peptide)
    print(lipdata)

print(timepoint_2_num_sig_peptides)
print(LipMS_all_data)
if len(u_timepoints) < 2:
    print(f'NOT ENOUGH SIG LiPMS peptides at multiple time points. Only timepoints with reported sig peptides are {u_timepoints}. Need atleast 2 unique time points, Exitting...')
    quit()
#load rep_ent data
rep_ent_file = sys.argv[1]
with open(rep_ent_file, 'rb') as fh:
    state_ent_data = pickle.load(fh)
print(list(state_ent_data.keys()))
#['u_ent_states_keys', 'u_ent_states_rep_ents', 'u_ent_states_frames']
pvalues = []
global obs_overlap_consistency, obs_exposure_consistency, obs_has_sig_overlap
for state_id, state_key in state_ent_data['u_ent_states_keys'].items():
    print('\n-----------------------------------------------------')
    rep_ent_data = load_rep_ent_data(state_ent_data['u_ent_states_rep_ents'][state_id])

    obs_overlap_consistency, obs_exposure_consistency = consistency(LipMS_sig_data, u_timepoints, rep_ent_data, 'O')
    print(obs_overlap_consistency, obs_exposure_consistency)

    obs_has_sig_overlap = False
    found_overlap = 0
    found_overlap_at_max = False
    for timepoint_idx, timepoint in enumerate(u_timepoints):
        if obs_overlap_consistency[timepoint_idx] > 0:
            found_overlap += 1
            if timepoint == max_timepoint:
                found_overlap_at_max = True
    if found_overlap >= 2 and found_overlap_at_max == True:
        obs_has_sig_overlap = True
    else:
        obs_has_sig_overlap = False
    pvalue = 0

    pvalue_counts = list(Parallel(n_jobs=nproc)(delayed(consistency)(LipMS_all_data, u_timepoints, rep_ent_data, 'T') for it  in range(reps)))
    pvalue = np.mean(pvalue_counts)
    print(f'pvalue: {pvalue} | state: {state_id}')
    pvalues += [pvalue]

rejected, q_value = fdrcorrection(pvalues)
print(rejected)
print(q_value)

state_ent_data['pvalues'] = {}
for state_id, state_key in state_ent_data['u_ent_states_keys'].items():
    state_ent_data['pvalues'][state_id] = q_value[state_id]

with open(f'{outfile}_consistency.pkl', 'wb') as fh:
    pickle.dump(state_ent_data, fh)
print(f'SAVED: {outfile}_consistency.pkl')

#END
end_time = time.time() - start_time
print(f'NORAML TERMINATION: {end_time}')
