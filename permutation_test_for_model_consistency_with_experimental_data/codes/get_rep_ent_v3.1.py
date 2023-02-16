#!/usr/bin/env python3

import numpy as np
import sys,os
import matplotlib
import glob
import pickle
import pandas as pd
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

np.set_printoptions(precision=3)

if len(sys.argv) != 6:
    print('[1] path to pkl file')
    print('[2] start frame')
    print('[3] end frame')
    print('[4] outfile basename')
    print('[5] path to secondary structure file')
    sys.exit()


def load_data(file_path):
    #load data
    for frame in range(start, end+1):
        #print(frame)

        for i,(k,v) in enumerate(orig_data[frame].copy().items()):
            crossings = []
            surr = []

            if frame != -1:
                #print(i,k,v)

                #Nterminal changbe
                if len(v[1][0]) > 0:

                    change_type = v[1][0][0]
                    ref_g = v[1][0][1]
                    traj_g = v[1][0][2]
                    #print(change_type, ref_g, traj_g)

                    if change_type in [0,1]: #gain
                        if np.abs(traj_g) < 0.6:
                            orig_data[frame][k][1][0] = np.asarray([])
                        else:
                            orig_data[frame][k][1][0] = np.asarray(orig_data[frame][k][1][0])
                            crossings += v[2]
                            surr += v[-1]

                    elif change_type in [2,3]: #loss
                        if np.abs(ref_g) < 0.6:
                            orig_data[frame][k][1][0] = np.asarray([])
                        else:
                            orig_data[frame][k][1][0] = np.asarray(orig_data[frame][k][1][0])
                            crossings += v[2]
                            surr += v[-1]

                    elif change_type in [4]: #switch chirl
                        if np.abs(traj_g) < 0.6 or np.abs(ref_g) < 0.6:
                            orig_data[frame][k][1][0] = np.asarray([])
                        else:
                            orig_data[frame][k][1][0] = np.asarray(orig_data[frame][k][1][0])
                            crossings += v[2]
                            surr += v[-1]
                else:
                    orig_data[frame][k][1][0] = np.asarray([])


                if len(v[1][1]) > 0:

                    change_type = v[1][1][0]
                    ref_g = v[1][1][1]
                    traj_g = v[1][1][2]
                    #print(change_type, ref_g, traj_g)

                    if change_type in [0,1]: #gain
                        if np.abs(traj_g) < 0.6:
                            orig_data[frame][k][1][1] = np.asarray([])
                        else:
                            orig_data[frame][k][1][1] = np.asarray(orig_data[frame][k][1][1])
                            v[1][1] = np.asarray(v[1][1])
                            crossings += v[2]
                            surr += v[-1]

                    if change_type in [2,3]: #loss
                        if np.abs(ref_g) < 0.6:
                            orig_data[frame][k][1][1] = np.asarray([])
                        else:
                            orig_data[frame][k][1][1] = np.asarray(orig_data[frame][k][1][1])
                            crossings += v[2]
                            surr += v[-1]

                    if change_type in [4]: #switch chirl
                        if np.abs(traj_g) < 0.6 or np.abs(ref_g) < 0.6:
                            orig_data[frame][k][1][1] = np.asarray([])
                        else:
                            orig_data[frame][k][1][1] = np.asarray(orig_data[frame][k][1][1])
                            crossings += v[2]
                            surr += v[-1]

                else:
                    orig_data[frame][k][1][1] = np.asarray([])

            #check if this native contact should be ignored
            if (orig_data[frame][k][1][0].size == 0 and orig_data[frame][k][1][1].size == 0) or (len(crossings) == 0):
                del orig_data[frame][k]

    return orig_data


def cluster(data, sec_elems):

    cluster_dict = {}
    cluster_rep_dict = {}

    frames = np.arange(start, end+1)
    print(f'frames to cluster: {frames}')


    #cluster based on (num_Ncross, num_Ccross, gN, gC, changeN, changeC)
    for frame in frames:
        if frame != -1 and frame in data:

            ent_data = data[frame]
            #print(frame, ent_data)
            #if frame not in cluster_dict:
            #    cluster_dict[frame] = {}

            for nc, ent in ent_data.items():
                print('\n', frame, nc, ent)
                #84, 123) [array([0.73461379, 0.11097391]), [[0, -0.3360524409177605, 0.7346137865365975], []], [[44]], [[41, 42, 43, 45, 46, 90, 91, 92, 113, 114, 115, 117, 118, 121]]]

                Ncross = tuple(ent[2][0])
                Ccross = tuple(ent[2][1])
                num_Ncross = int(len(Ncross))
                num_Ccross = int(len(Ccross))

                if num_Ncross == 0 and num_Ccross == 0:
                    continue

                if len(ent[1][0]) > 0:
                    gN = np.round(ent[0][0]).astype(int)
                    changeN = int(ent[1][0][0])

                else:
                    gN = 99
                    changeN = -1

                if len(ent[1][1]) > 0:
                    gC = np.round(ent[0][1]).astype(int)
                    changeC = int(ent[1][1][0])

                else:
                    gC = 99
                    changeC = -1

                #key = (firstNcross_2strcut, firstCcross_2strcut, num_Ncross, num_Ccross, gN, gC, changeN, changeC)
                #key = (N2structs_dist, C2structs_dist, num_Ncross, num_Ccross, gN, gC, changeN, changeC)
                key = (num_Ncross, num_Ccross, gN, gC, changeN, changeC, Ncross, Ccross)
                print(key)

                if key not in cluster_dict:
                    cluster_dict[key] = {}

                if frame not in cluster_dict[key]:
                    cluster_dict[key][frame] = {}

                cluster_dict[key][frame][nc] = ent


    #condense cluster dict by combining keys that have the same num_Ncross, num_Ccross, gN, gC, changeN, changeC and their Ncross and Ccross are within +/- 10 residues of eachother
    cluster_dict_copy = cluster_dict.copy()
    already_checked_keys = []
    for key1,frames1 in cluster_dict_copy.items():
        print(f'\nKEY1: {key1}')
        #print(frames1.keys())

        #check if first 6 parameters are the same
        for key2,frames2 in cluster_dict_copy.items():

            if key2 not in already_checked_keys:
                #print(f'KEY2: {key2}')
                #print(frames2.keys())

                if key1 != key2:
                    if key1[0:6] == key2[0:6]:
                        Npass = False
                        Cpass = False

                        #check N term crossings
                        if key1[6] != ():
                            key1_Nset = set(key1[6])
                        else:
                            key1_Nset = set()
                        #print(key1_Nset)

                        if key2[6] != ():
                            key2_Nset = set(np.hstack([np.arange(x-10, x+11) for x in key2[6]]))
                        else:
                            key2_Nset = set()
                        #print(key2_Nset)

                        if key1[7] != ():
                            key1_Cset = set(key1[7])
                        else:
                            key1_Cset = set()
                        #print(key1_Cset)

                        if key2[7] != ():
                            key2_Cset = set(np.hstack([np.arange(x-10, x+11) for x in key2[7]]))
                        else:
                            key2_Cset = set()
                        #print(key2_Cset)


                        if key1_Nset == key2_Nset:
                            Npass = True
                        if key1_Cset == key2_Cset:
                            Cpass = True

                        if len(key1_Nset.intersection(key2_Nset)) > 0:
                            Npass = True
                        if len(key1_Cset.intersection(key2_Cset)) > 0:
                            Cpass = True

                        if Npass == True and Cpass == True:
                            #print(Npass, Cpass)

                            for frame,nc in cluster_dict[key1].items():
                                if frame in cluster_dict[key2]:
                                    cluster_dict[key2][frame].update(cluster_dict[key1][frame])
                                else:
                                    cluster_dict[key2][frame] = cluster_dict[key1][frame]
                            del cluster_dict[key1]
                            already_checked_keys += [key1]
                            break
                else:
                    continue


    #########################################################################
    # make representative entanglement for each cluster in each frame
    print('# make representative entanglement for each cluster in each frame')
    rep_dict = {}
    #print(cluster_dict.keys())

    frame2ent_dict = {}
    for key,frames in cluster_dict.items():
        print(f'\nKEY: {key}')
        #print(frames.keys())

        rep_dict[key] = {}

        min_loop = 99999
        for frame,ents in frames.items():
            print(f'FRAME: {frame}')

            if frame not in frame2ent_dict:
                frame2ent_dict[frame] = [key]
            else:
                frame2ent_dict[frame] += [key]

            for nc, ent in ents.items():
                print(nc, ent)

                loop_l = np.diff(nc)[0]
                #print(loop_l)

                if loop_l < min_loop:
                    min_loop = loop_l
                    rep_nc = nc
                    rep_ent = ent
                    rep_frame = frame

        print('REP:', key, rep_nc, [rep_frame, rep_ent])
        rep_dict[key][rep_nc] = [rep_frame, rep_ent]

    #########################################################################
    #for each frame determine the unique ent present to determine the ent_state of that frame
    for frame in sorted(frame2ent_dict.keys()):
        keys = frame2ent_dict[frame]
        print(keys)
        print(set(keys))
        #keys = np.vstack(keys).astype(float)
        frame2ent_dict[frame] = set(keys)
        print(frame)
        print(frame2ent_dict[frame])

    # for each unique ent_state across all frames get the rep_ent for each key in state and all frames associated with state
    print('\nfor each unique ent_state across all frames get the rep_ent for each key in state and all frames associated with state')
    u_ent_states = {}
    u_ent_states_rep_ents = {}
    u_ent_states_frames = {}
    u_ent_states_master = {}

    #print(frame2ent_dict.values(), len(list(frame2ent_dict.values())))
    #u_states = {array.tostring(): array for array in frame2ent_dict.values()}
    #u_states = u_states.values()
    u_states = np.unique(list(frame2ent_dict.values()))
    #print(u_states, len(u_states))
    for u_idx, u_ent_state in enumerate(u_states):

        u_ent_state_rep_ent = []
        u_ent_state_frames = []
        for key in u_ent_state:
            rep_ent = rep_dict[tuple(key)]
            u_ent_state_rep_ent += [rep_ent]

        for frame,ent_state in frame2ent_dict.items():
            if np.array_equal(ent_state, u_ent_state):
                u_ent_state_frames += [frame]

        u_ent_states[u_idx] = u_ent_state
        u_ent_states_rep_ents[u_idx] = u_ent_state_rep_ent
        u_ent_states_frames[u_idx] = u_ent_state_frames
        print(f'\nUnique Entangled State: {u_idx}')
        print(u_ent_state)
        for u_rep in u_ent_state_rep_ent:
            print(u_rep)
        print(u_ent_state_frames)

    #then save everything to a master dictionary
    u_ent_states_master['u_ent_states_keys'] = u_ent_states
    u_ent_states_master['u_ent_states_rep_ents'] = u_ent_states_rep_ents
    u_ent_states_master['u_ent_states_frames'] = u_ent_states_frames

    return cluster_dict, rep_dict, u_ent_states_master

#load user arguments
global tart, end, aID
start = int(sys.argv[2])
end = int(sys.argv[3])
outfile_basename = sys.argv[4]
start_time = time.time()

files = glob.glob(sys.argv[1])[0]

#make a matrix containing the secondary structure elements
sec_elems_file_path = sys.argv[5]
if sec_elems_file_path != 'nan':
    sec_elems = np.asarray([[int(row[1]), int(row[2])] for row in [x.strip('\n').split() for x  in open(sec_elems_file_path).readlines()]])
    sec_elems = np.asarray(sorted(sec_elems, key=lambda x:x[0]))
    print(f'sec_elems_file_path: {sec_elems_file_path}')
    print(f'sec_elems: {sec_elems}')

with open(files,'rb') as fh:
    orig_data = pickle.load(fh)

#check if end supplied is greater than the number of frames
if end > len(orig_data):
    print(f'frames in {files}')
    print(orig_data.keys())
    num_frames = len([x for x in orig_data.keys() if x != -1])
    end_frame = max(orig_data.keys())
    print(f'supplied end frame {end} > amount of frames in .pkl file {num_frames}\nsetting end={end_frame}')
    end = end_frame

#load in orgiinal data
orig_data = load_data(orig_data)

#cluster changes in ent
cluster_data, rep_data, u_ent_states_data = cluster(orig_data, sec_elems)
with open(f'{outfile_basename}_raw_cluster.pkl', 'wb') as fh:
    pickle.dump(cluster_data, fh)
print(f'SAVED: {outfile_basename}_rep_ent_raw_cluster.pkl')

with open(f'{outfile_basename}_u_ent_states.pkl', 'wb') as fh:
    pickle.dump(u_ent_states_data, fh)
print(f'SAVED: {outfile_basename}_u_ent_states.pkl')

print('--------------------------------------------------------------------------------------------------------------')
print('\nRepresentative changes in ent from clustering')
outdict = {'rep_ent_idx': [], 'frame': [], 'NC': [], 'Nchange': [], 'Nchange_ref_g': [], 'Nchange_frame_g': [], 'Cchange': [], 'Cchange_ref_g': [], 'Cchange_frame_g': [], 'crossing': [], 'surrounding': []}
for rep_i, (key, ent_info) in enumerate(rep_data.items()):
    #(#cross, Nchir, Cchir, Nchange, Cchange)
    print(f'\nRep Ent Number: {rep_i}')
    #print(f'Nearest 2struct Ncrossings: {key[0]}')
    #print(f'Nearest 2struct Ccrossings: {key[1]}')
    print(f'N Number of Crossings: {key[0]}')
    print(f'C Number of Crossings: {key[1]}')
    print(f'gN: {key[2]}')
    print(f'gC: {key[3]}')
    print(f'N ent change type: {key[4]}')
    print(f'C ent change type: {key[5]}')
    for k,(frame, ((gNvalue, gCvalue), (gNchange, gCchange), crossings, surr)) in ent_info.items():
        print(k, frame, gNvalue, gCvalue, gNchange, gCchange , crossings, surr)

        print(f'Frame: {frame}')
        print(f'Native contact: {k}')
        outdict['rep_ent_idx'] += [rep_i]
        outdict['frame'] += [frame]
        outdict['NC'] += [k]

        if gNchange.size > 0:
            print(f'N-term frame partial linking value: {gNvalue}')
            print(f'N-term ref partial linking value: {gNchange[1]}')
            outdict['Nchange'] += [key[4]]
            outdict['Nchange_ref_g'] += [gNchange[1]]
            outdict['Nchange_frame_g'] += [gNvalue]
        else:
            outdict['Nchange'] += [np.nan]
            outdict['Nchange_ref_g'] += [np.nan]
            outdict['Nchange_frame_g'] += [np.nan]

        if gCchange.size > 0:
            print(f'C-term frame partial linking value: {gCvalue}')
            print(f'C-term ref partial linking value: {gCchange[1]}')
            outdict['Cchange'] += [key[5]]
            outdict['Cchange_ref_g'] += [gCchange[1]]
            outdict['Cchange_frame_g'] += [gCvalue]
        else:
            outdict['Cchange'] += [np.nan]
            outdict['Cchange_ref_g'] += [np.nan]
            outdict['Cchange_frame_g'] += [np.nan]

        print(f'Crossings: {crossings}')
        print(f'Residues with 8A of crossings: {surr}')

        outdict['crossing'] += [crossings]
        outdict['surrounding'] += [surr]

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1500)
df = pd.DataFrame(data=outdict)
print(df)
df.to_csv(f'{outfile_basename}_summary_rep_ents.csv', float_format='%.3f')
print(f'SAVED: {outfile_basename}_summary_rep_ents.csv')

print(f'NORAML TERMINATION @ {time.time() - start_time}')
