#!/usr/bin/env python3
################################################################################################################
script_title='entanglement_analysis'
script_version=2.0
script_author='Ian Sitarik'
script_updatelog=f"""Update log for {script_title} version {script_version}

                   Date: 08.20.2021
                   Note: Started covertion of topology_anal codes

                   Date: 09.01.2021
                   Note: v1.5 reproduced paper. This version will include a separate feature to iterate through all loops

                   Date: 09.11.2021
                   note: implimented joblib for multiprocessing support ~O(|native contacts|)

                   Date: 10.06.2021
                   note: added ability to track resid from PDB for entanglement

                   Date: 10.24.2021
                   note:- calculates two change in entanglement methods:
                        1. if there is a change in the un rounded partial linking number of any tail for a given
                        NC greater than the threshold provided by the user
                        2. if the total linking number changes

                        - added the ability to catch missing residues in two ways.
                        1. if there are any missing residues between the first and last residue in the PDB
                        2. from the PDB file its self as they should be reported unless the user removed the info
                            as happens sometimes

                        NOTE::: make sure reference and traj have same resid labeling

                   Date: 10.28.2021
                   note: added automatic restart functionaily

                   Date: 10.30.2021
                   note: added control file read in

                   Date: 10.30.2021
                   note: added fraction of native contacts analysis (automatic)

                   Date: 11.17.2021
                   note: added minimal loop finder

                   Date: 11.22.2021
                   note: added mask to select certain part of PDB

                   Date: 11.24.2021
                   note: added correction to output loop_analysis for ref state as well
                   note: added unique entanglement finder

                   Date: 02.08.2022
                   note: REBIRTH of script from OP_v8.0.py
                   note: overhaul to simplify inputs and calculations

                   Date: 03.05.2022
                   note: simplified output to two files
                   1. the time series for fraction of native contacts (Q) and the fraction of native contacts with a change in entanglement (G)
                   2. a pickle binary file containing a single dictionary with one entry for the reference state and a single entry for each frame analyzed in the trajectory
                        this file can be examined by launching an interactive python session and using the following commands

                        import pickle
                        with open('./test_outfiles/entanglement_analysis_1.4/output/6u32.pkl', 'rb') as fh:
                            data = pickle.load(fh)

                        for k,v in data.items():
                            print(k,v)


                        the top level of keys will be integers ranging from -1 to the number of frames analyzed minus 1. For eample if you
                        analyzed a trajectory with 10 frames the dictionary would have a total of 11 entries with the following keys

                        -1 = reference state results
                        0 = first frame results
                        1 = second frame results
                        ...
                        9 = tenth frame results


                        in each of these entires the value is another dictionary containing one entry for each native contact that was detected to have a change in entanglement

                        for the refernce state key = -1 the inner dictionary will be structured like this
                        key = (residues invovled in native contact)
                        value = [array containing partial linking values for the N and C termini for this native contact,
                                 array containing the N and C terminal crossings for this native contact,
                                 residues within 8A of the crossing residues]

                        so for example if the key value pair in the reference state returned this
                        (4, 50) [array([0.        , 0.84160559]), [[], [61]], [[], [24, 25, 26, 27, 28, 29, 34, 35, 36, 37, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63]]]
                        the native contact formed by the residues 4 and 50 has an entanglement present.
                        the partial linking value of the N terminus is 0 indicating no entanglement
                        the partial linking value of the C terminus is 0.84160559 indicating a entanglement is present
                        Therefor the N terminus should not have a crossing while the C terminus should and indeed does at residue 61
                        The residues who have alpha carbons within 8A of the crossing residues are reported last


                        for the frames in the traj key >= 0 the inner dictionary will be structured like this
                        key = (residues invovled in native contact)
                        value = [array containing partial linking values for the N and C termini for this native contact,
                                 array containing [change type, refernce state partial linking value, frame partial linking value] for the N and C termini for this native contact,
                                 array containing the N and C terminal crossings for this native contact,
                                 residues within 8A of the crossing residues]

                                change_type = 0 if gain of entanglement and no change in chirality
                                change_type = 1 if gain of entanglement and change in chirality
                                change_type = 2 if loss of entanglement and no change in chirality
                                change_type = 3 if loss of entanglement and change in chirality
                                change_type = 4 if only chirality

                        so for example if the key value pair in a frame returned this
                        (108, 135) [array([0.96495744, 0.00888425]), [[0, 0.06142746756886508, 0.9649574387268217], []], [[12, 106]], [[10, 11, 13, 14, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 131, 132, 133, 134, 135, 163, 164, 167]]]
                        the native contact formed by the residues 108 and 135 has an entanglement present.
                        the partial linking value of the N terminus is 0.96495744 indicating a entanglement
                        the partial linking value of the C terminus is 0.00888425 indicating no entanglement is present
                        This signifies a gain of entanglement with no change in chirality so the change_type = 0
                        and the patial linking value in the reference state for the N termini is 0.06142746756886508 while in the misfolded state is 0.9649574387268217
                        Therefor the C terminus should not have a crossing while the N terminus should and indeed does at residue 12 and 106
                        The residues who have alpha carbons within 8A of the crossing residues are reported last


                   Date: 06.15.2022
                   note: added ability to detect missing residues in reference state files and correct the entanglement output to remove entanglements if
                         1. there are any missing residues within +/- 10 of any crossing
                         2. more than three consecutive missing residues in the loop
                         3. more than 5% of the loop has missing residues

                         This should not affect the trajectory caclulations as those structure should be complete and a complete reference state should be used when analyzing them.
                         This only affects the case when you are analyzing a single structure for entanglements and not comparing to a treatment structure or trajectory.

                    Date: 08.16.2022
                    note: changed from finding surr residues using mdanalysis to manual calculation

                    Date: 08.19.2022
                    note: fixed error in fraction of native contact analysis that did not include 1.2 cutoff for thermal flux

                    Date: 11.21.2022
                    note: added ability to skip topoly if you do not need crossings or surrounding residues and only need GQ output

                  """

################################################################################################################

import os
import sys
import numpy as np
import time
import itertools
from MDAnalysis import *
from scipy.spatial.distance import pdist, squareform
from itertools import product, combinations
from itertools import chain as iterchain
from joblib import Parallel, delayed
import configparser
import pickle
from topoly import *
import more_itertools as mit

##################################################################################################################
### START argument parse ###
##################################################################################################################
if len(sys.argv) != 12:
    print(f'[1] = path to control file')
    print(f'[2] = outfile_basename')
    print(f'[3] = start_frame')
    print(f'[4] = end_frame')
    print(f'[5] = frame_stride')
    print(f'[6] = num processors')
    print(f'[7] = ref coor file')
    print(f'[8] = dcd file')
    print(f'[9] = psf')
    print(f'[10] = path to sec_structs file')
    print(f'[11] = use topoly False=No True=Yes')
    quit()


###### CONTROL FILE example ######
# ref_path = path to reference coordinate file (PDB or CRD)
# psf = path to protein structure file
# out_path = path to output directory
# sec_elems_file = path to sec elements file or nan if you want to skip Q
# #################################

config = configparser.ConfigParser()
config.read(sys.argv[1])
default = config['DEFAULT']
print(default)

in_path=sys.argv[8]
if in_path == None:
    print(f'\n{script_updatelog}\n')
    sys.exit()

elif in_path == 'nan':
    print(f'in_path: {in_path}')

else:
    in_paths=[x.strip('\n').split(', ') for x in os.popen(f'ls -v {in_path}').readlines()]
    for ipath in in_paths:
        print(f'{ipath}')

out_path=default['out_path']
print(f'out_path: {out_path}')
if out_path == None: print(f'\n{script_updatelog}\n'); sys.exit()

psf=sys.argv[9]
print(f'psf: {psf}')
if psf == None: print(f'\n{script_updatelog}\n'); sys.exit()

global use_topoly
use_topoly=sys.argv[11]
print(f'use_topoly: {use_topoly}')
if use_topoly == None: print(f'\n{script_updatelog}\n'); sys.exit()
if use_topoly == 'True': use_topoly=True
if use_topoly == 'False': use_topoly=False


print_updatelog=True
print(f'print_updatelog: {print_updatelog}')
if print_updatelog != None: print(f'\n{script_updatelog}\n')

global print_framesummary
print_framesummary=False
print(f'print_framesummary: {print_framesummary}')
if print_framesummary == None: print(f'\n{script_updatelog}\n')

outfile_basename=sys.argv[2]
print(f'outfile_basename: {outfile_basename}')
if outfile_basename == None: print(f'\n{script_updatelog}\n'); sys.exit()

start_frame=sys.argv[3]
print(f'start_frame: {start_frame}')
if start_frame == None:
    print(f'\n{script_updatelog}\n')
    start_frame = 0
else: start_frame=int(start_frame)

end_frame=sys.argv[4]
print(f'end_frame: {end_frame}')
if end_frame == None:
    print(f'\n{script_updatelog}\n')
    end_frame = 999999999999
else: end_frame=int(end_frame)

frame_stride=sys.argv[5]
print(f'frame_stride: {frame_stride}')
if frame_stride == None:
    print(f'\n{script_updatelog}\n')
    frame_stride = 1
else: frame_stride=int(frame_stride)

global S
S = 6
print(f'S: {S}')

global nproc
nproc = int(sys.argv[6])
print(f'nproc: {nproc}')

global change_threshold
change_threshold = 0.0
print(f'change_threshold: {change_threshold}')

global ent_gval_threshold
ent_gval_threshold = 0.6
print(f'ent_gval_threshold: {ent_gval_threshold}')

ref_path= sys.argv[7]
print(f'ref_path: {ref_path}')
if ref_path == None: print(f'\n{script_updatelog}\n'); sys.exit()

ref_mask = default['ref_mask']
print(f'ref_mask: {ref_mask}')
if ref_mask == None: ref_mask = 'all'

traj_mask = default['traj_mask']
print(f'traj_mask: {traj_mask}')
if traj_mask == None: traj_mask = 'all'

global TS
TS = 5
print(f'TS: {TS}')

#make a matrix containing the secondary structure elements
sec_elems_file_path = sys.argv[10]
if sec_elems_file_path != 'nan':
    sec_elems = np.asarray([[int(row[1]), int(row[2])] for row in [x.strip('\n').split() for x  in open(sec_elems_file_path).readlines()]])
    print(f'sec_elems_file_path: {sec_elems_file_path}')
    print(f'sec_elems: {sec_elems}')

##################################################################################################################
### START initial loading of structure files and qualtiy control ###
##################################################################################################################

### START dir declaration ###

if os.path.exists(f'{out_path}/'):
    print(f'{out_path}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/'):
    print(f'{out_path}{script_title}_{script_version}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/logs/'):
    print(f'{out_path}{script_title}_{script_version}/logs/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/logs/')

if os.path.exists(f'{out_path}{script_title}_{script_version}/output/'):
    print(f'{out_path}{script_title}_{script_version}/output/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}_{script_version}/output/')

### END dir declaration ###

### START preference setting ###

start_time=time.time() #time since epoch
print('time since epoch = '+str(start_time))

np.set_printoptions(precision=4, suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)
np.seterr(divide='ignore')

### END preference setting ###

######################################################################################################################
# USER DEFINED FUNCTIONS                                                                                             #
######################################################################################################################

def gen_nc_gdict(coor, coor_cmap, **kwargs):
    dom_nc_gdict = {}
    dom_gn_dict = {}
    dom_contact_ent = {}
    global dot_matrix, l

    nc_indexs = np.stack(np.nonzero(coor_cmap)).transpose()

    l = len(coor)
    #print(f'l: {l}')

    #make R and dR waves of length N-1
    range_l = np.arange(0, l-1)
    range_next_l = np.arange(1,l)

    coor = coor.astype(np.float32)
    R = 0.5*(coor[range_l] + coor[range_next_l])
    dR = coor[range_next_l] - coor[range_l]

    #make dRcross matrix
    pair_array = np.asarray(list(itertools.product(dR,dR)))

    x = pair_array[:,0,:]
    y = pair_array[:,1,:]

    dR_cross = np.cross(x,y)

    #make Rnorm matrix
    pair_array = np.asarray(list(itertools.product(R,R)))

    diff = pair_array[:,0,:] - pair_array[:,1,:]
    diff = diff.astype(np.float32)
    Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3
    Runit = Runit.astype(np.float32)

    #make final dot matrix
    dot_matrix = [np.dot(x,y) for x,y in zip(Runit,dR_cross)]
    dot_matrix = np.asarray(dot_matrix)
    dot_matrix = dot_matrix.reshape((l-1,l-1))

    coor = [list(x) for x in coor]
    contact_ent = list(Parallel(n_jobs=nproc)(delayed(g_calc)(i, j, coor) for i,j in nc_indexs if j >= i + 10))
    contact_ent = {k: v for d in contact_ent for k, v in d.items()}

    return (contact_ent)


def g_calc(i,j, coor):

    loop_range = np.arange(i,j)
    crossings = [list(), list()]

    nterm_range = np.arange(TS,i-S)
    gn_parital_link, gn, gn_j1, gn_j2 = sub_lists(nterm_range,loop_range)
    if abs(gn_parital_link) > ent_gval_threshold:

        if use_topoly == True:
            Ncrossings = lasso_type(coor, [(i,j)], more_info=True)
            #crossings[0] = Ncrossings[(i,j)]['crossingsN']
            crossings[0] = list(np.abs(np.unique(np.asarray([x.replace('*','') for x in Ncrossings[(i,j)]['crossingsN']], dtype=int))))
        elif use_topoly == False:
            pass

    cterm_range = np.arange(j+S,l-TS-1)
    gc_parital_link, gc, gc_j1, gc_j2 = sub_lists(cterm_range,loop_range)
    if abs(gc_parital_link) > ent_gval_threshold:

        if use_topoly == True:
            Ccrossings = lasso_type(coor, [(i,j)], more_info=True)
            crossings[1] = list(np.abs(np.unique(np.asarray([x.replace('*','') for x in Ccrossings[(i,j)]['crossingsC']], dtype=int))))
            #crossings[1] = Ccrossings[(i,j)]['crossingsC']

        elif use_topoly == False:
            pass

    #print(i,j,gn_parital_link, gc_parital_link, crossings)
    out = {(i, j):np.asarray([[gn_parital_link,gc_parital_link],crossings])}

    return out



#@nb.njit(fastmath=True)
def helper_func(g_vals: np.ndarray):

    return abs(g_vals.sum()/(4.0*np.pi))


def sub_lists(thread, loop):

    if len(thread) > 10:
        pairs = np.fromiter(itertools.chain(*itertools.product(thread, loop)), int).reshape(-1, 2)
        parital_link = dot_matrix[pairs[:,0], pairs[:,1]].sum()/(4.0*np.pi)

        return parital_link, parital_link, thread[0], thread[-1]

    else:
        return 0, 0, 0, 0


def ent_cmap(cor, ref = True, restricted = True, cut_off = 8.0, bb_buffer = 4, **kwargs):
    if print_framesummary: print(f'\nCMAP generator')
    if print_framesummary: print(f'ref: {ref}\nrestricted: {restricted}\ncut_off: {cut_off}\nbb_buffer: {bb_buffer}')

    distance_map=squareform(pdist(cor,'euclidean'))
    distance_map=np.triu(distance_map,k=bb_buffer)

    contact_map = np.where((distance_map < cut_off) & (distance_map > 0), 1, 0)

    contact_num=contact_map.sum()

    if print_framesummary: print(f'Total number of contacts in is {contact_num}')

    return contact_map, contact_num

def change_anal(frame_gln_data, ref_gln_data):
    Nterm_chng_ent_dict={}
    Cterm_chng_ent_dict={}
    G_counter = 0

    for nc,ent_info in ref_gln_data.items():
        gvals, _, _ = ent_info

        if nc in frame_gln_data:

            frame_gvals, _, _ = frame_gln_data[nc]
            change_found = 0
            for tail_idx in [0,1]:

                frame_g = frame_gvals[tail_idx]
                ref_g = gvals[tail_idx]
                if abs(frame_g) > ent_gval_threshold or abs(ref_g) > ent_gval_threshold:
                    if abs(frame_g - ref_g) > change_threshold:
                        if abs(frame_g.round()) > abs(ref_g.round()) and frame_g.round()*ref_g.round() >= 0:
                            change_found = 1
                            #print('gain no chirality shift')
                            if tail_idx == 0:
                                Nterm_chng_ent_dict[nc] = [0, ref_g, frame_g]
                            if tail_idx == 1:
                                Cterm_chng_ent_dict[nc] = [0, ref_g, frame_g]

                        if abs(frame_g.round()) > abs(ref_g.round()) and frame_g.round()*ref_g.round() < 0:
                            change_found = 1
                            #print('gain and chirality shift')
                            if tail_idx == 0:
                                Nterm_chng_ent_dict[nc] = [1, ref_g, frame_g]
                            if tail_idx == 1:
                                Cterm_chng_ent_dict[nc] = [1, ref_g, frame_g]

                        if abs(frame_g.round()) < abs(ref_g.round()) and frame_g.round()*ref_g.round() >= 0:
                            change_found = 1
                            #print('loss no chirality shift')
                            if tail_idx == 0:
                                Nterm_chng_ent_dict[nc] = [2, ref_g, frame_g]
                            if tail_idx == 1:
                                Cterm_chng_ent_dict[nc] = [2, ref_g, frame_g]

                        if abs(frame_g.round()) < abs(ref_g.round()) and frame_g.round()*ref_g.round() < 0:
                            change_found = 1
                            #print('loss and chirality shift')
                            if tail_idx == 0:
                                Nterm_chng_ent_dict[nc] = [3, ref_g, frame_g]
                            if tail_idx == 1:
                                Cterm_chng_ent_dict[nc] = [3, ref_g, frame_g]

                        if abs(frame_g.round()) == abs(ref_g.round()) and frame_g*ref_g < 0 and frame_g.round() != 0 and ref_g.round() != 0:
                            change_found = 1
                            #print('pure chirality shift')
                            if tail_idx == 0:
                                Nterm_chng_ent_dict[nc] = [4, ref_g, frame_g]
                            if tail_idx == 1:
                                Cterm_chng_ent_dict[nc] = [4, ref_g, frame_g]

                    else:
                        continue
                else:
                    continue

            if change_found == 1:
                G_counter += 1

        else:
            continue


    return Nterm_chng_ent_dict, Cterm_chng_ent_dict, G_counter


######################################################################################################################
# MAIN                                                                                                               #
######################################################################################################################

global framenum,frametime
framenum = -1
frametime = 0.0
outdata = {}
outdata[framenum] ={}
#load ref structs and get entanglement
if ref_path != 'nan':
    print('\n########################################START loading of reference universe########################################\n')

    #load data into a universe and use the mask the suer provided
    print(f'Loading: {ref_path}')
    ref = Universe(ref_path)
    ref_calphas = ref.select_atoms(f'{ref_mask}')
    print(ref_calphas, len(ref_calphas))

    #get mapping for coor_idx to PDB resid
    global ref_cooridx2pdbresid
    ref_cooridx2pdbresid = {i:res.resid for i,res in enumerate(ref_calphas.residues)}
    ref_pdbresid2cooridx = {v: k for k, v in ref_cooridx2pdbresid.items()}
    #print(ref_cooridx2pdbresid)

    #get coordinate positions
    ref_coor = ref_calphas.positions
    ref_distance_map=squareform(pdist(ref_coor,'euclidean'))
    #print(ref_coor[:10])


    #get cmap for G and restricted cmap for Q if nec_elems file was specified
    ref_cmap, ref_num_contacts = ent_cmap(ref_coor)

    if sec_elems_file_path != 'nan':
        ref_distance_map=squareform(pdist(ref_coor,'euclidean'))
        ref_contact_map = np.zeros((ref_distance_map.shape))
        ref_contact_map[ref_distance_map < 8.0] = 1
        ref_contact_map=np.triu(ref_contact_map,k=4)

        restricted_ref_cmap = ref_contact_map.copy()

        sec_residues = []
        for elems in sec_elems:
            sec_residues.append(np.arange(elems[0]-1, elems[1]))
        sec_residues = np.hstack((sec_residues))
        anti_sec_residues = [x for x in np.arange(0,len(ref_coor)) if x not in sec_residues]

        restricted_ref_cmap[:,anti_sec_residues] = 0
        restricted_ref_cmap[anti_sec_residues,:] = 0
        restricted_ref_num_contacts = restricted_ref_cmap.sum()
        print(f'restricted_ref_num_contacts: {restricted_ref_num_contacts}')


    #GLN analysis
    ref_cont_ent_data = gen_nc_gdict(ref_coor, ref_cmap)
    for nc,(gvals, crossings) in ref_cont_ent_data.copy().items():
        #print(nc,gvals, crossings)

        del ref_cont_ent_data[nc]

        nc = tuple(ref_cooridx2pdbresid[n] for n in nc)
        Ncrossings = crossings[0]
        if len(Ncrossings) != 0:
            distance_map = np.where(ref_distance_map < 8.0, 1, 0)[Ncrossings, :]
            Nsurr = list(np.where(distance_map == 1)[1])
            Nsurr = [ref_cooridx2pdbresid[s] for s in Nsurr]
        else:
            Nsurr = []

        Ccrossings = crossings[1]
        if len(Ccrossings) != 0:
            distance_map = np.where(ref_distance_map < 8.0, 1, 0)[Ccrossings, :]
            Csurr = list(np.where(distance_map == 1)[1])
            Csurr = [ref_cooridx2pdbresid[s] for s in Csurr]
        else:
            Csurr = []

        Ncrossings = [ref_cooridx2pdbresid[c] for c in Ncrossings]
        Ccrossings = [ref_cooridx2pdbresid[c] for c in Ccrossings]
        #print(nc,gvals, crossings)
        ref_cont_ent_data[nc] = [gvals, [Ncrossings, Ccrossings], [Nsurr, Csurr]]
        outdata[framenum][nc] = [gvals, [Ncrossings, Ccrossings], [Nsurr, Csurr]]
        #print(nc, outdata[framenum][nc])

### END loading of analysis universe ###

### START analysis of universe ###
trajdata = {}
traj_crossings = {}
frame_times = []
GQ_list = []
print('\n########################################START analysis of trajectory########################################\n')

if in_path != 'nan':

    #get alpha carbons atoms and then positions of them
    print(f'Loading: {psf} & {in_paths}')
    u = Universe(psf,in_paths)
    #u = Universe(in_paths)
    u_calphas = u.select_atoms(f'{traj_mask}')

    print(u_calphas,  len(u_calphas))

    total_frames = len(u.trajectory[start_frame:end_frame:frame_stride])
    print(f'Total_frames to analyze: {total_frames}')
    for frame_idx,ts in enumerate(u.trajectory[start_frame:end_frame:frame_stride]):

        print(f'\n\nFrame: {ts.frame}')

        framenum = ts.frame
        frametime = np.around(ts.time, decimals=4)
        frame_start_time=time.time() #time since epoch

        #get mapping for coor_idx to PDB resid
        global frame_cooridx2pdbresid
        frame_cooridx2pdbresid = {i:res.resid for i,res in enumerate(u_calphas.residues)}
        frame_pdbresid2cooridx = {v: k for k, v in frame_cooridx2pdbresid.items()}
        #print(frame_cooridx2pdbresid)

        frame_coor = u_calphas.positions
        frame_distance_map=squareform(pdist(frame_coor,'euclidean'))
        #print(frame_coor[:10])

        #get cmap for G and restricted cmap for Q if nec_elems file was specified
        frame_cmap, frame_num_contacts = ent_cmap(frame_coor)
        frame_cmap = frame_cmap*ref_cmap

        Q = 0
        G = 0
        if sec_elems_file_path != 'nan':
            frame_distance_map=squareform(pdist(frame_coor,'euclidean'))
            frame_contact_map = np.zeros((frame_distance_map.shape))
            frame_contact_map[frame_distance_map < ref_distance_map*1.2] = 1
            frame_contact_map=np.triu(frame_contact_map,k=4)
            restricted_frame_cmap = frame_contact_map.copy()

            sec_residues = []
            for elems in sec_elems:
                sec_residues.append(np.arange(elems[0]-1, elems[1]))
            sec_residues = np.hstack((sec_residues))
            anti_sec_residues = [x for x in np.arange(0,len(frame_coor)) if x not in sec_residues]

            restricted_frame_cmap[:,anti_sec_residues] = 0
            restricted_frame_cmap[anti_sec_residues,:] = 0
            restricted_frame_cmap = restricted_frame_cmap*restricted_ref_cmap
            restricted_frame_num_contacts = restricted_frame_cmap.sum()
            #print(f'restricted_frame_num_contacts: {restricted_frame_num_contacts}')
            Q = restricted_frame_num_contacts/restricted_ref_num_contacts
            print(f'Q: {Q}')

        #GLN analysis
        frame_cont_ent_data = gen_nc_gdict(frame_coor, frame_cmap)

        for nc,(gvals, crossings) in frame_cont_ent_data.copy().items():
            #print(nc,gvals, crossings)

            nc = tuple(frame_cooridx2pdbresid[n] for n in nc)
            Ncrossings = crossings[0]
            if len(Ncrossings) != 0:
                distance_map = np.where(frame_distance_map < 8.0, 1, 0)[Ncrossings, :]
                Nsurr = list(np.where(distance_map == 1)[1])
                Nsurr = [frame_cooridx2pdbresid[s] for s in Nsurr]
            else:
                Nsurr = []

            Ccrossings = crossings[1]
            if len(Ccrossings) != 0:
                distance_map = np.where(frame_distance_map < 8.0, 1, 0)[Ccrossings, :]
                Csurr = list(np.where(distance_map == 1)[1])
                Csurr = [frame_cooridx2pdbresid[s] for s in Csurr]
            else:
                Csurr = []

            Ncrossings = [frame_cooridx2pdbresid[c] for c in Ncrossings]
            Ccrossings = [frame_cooridx2pdbresid[c] for c in Ccrossings]
            #print(nc,gvals, crossings)
            frame_cont_ent_data[nc] = [gvals, [Ncrossings, Ccrossings], [Nsurr, Csurr]]
            #print(nc, frame_cont_ent_data[nc])


        #change in entanglement analysis
        #resulting dict(nc:[change_type, tail_idx, ref_g, frame_g])
        #where change_type 0:gain not chirality shift, 1:gain and chirality shift
        #where change_type 2:loss not chirality shift, 3:lost and chirality shift
        #where change_type 4:pure chirality shift
        #where tail_idx 0:Nterm 1:Cterm
        frame_Nterm_change_results, frame_Cterm_change_results, G = change_anal(frame_cont_ent_data, ref_cont_ent_data)

        #calculate number of native contacts with a change in entanglement from ref state
        change_ncs = list(frame_Nterm_change_results.keys()) + list(frame_Cterm_change_results.keys())
        G = G/frame_cmap.sum()
        print(f'G: {G}')
        GQ_list.append([framenum, frametime, G, Q])

        #generate outdata
        outdata[framenum] = {}
        #find residues within 8A of crossing


        #loop through native contacts that have a change in entanglement and add them to output dictionary
        for nc in change_ncs:

            if nc in outdata[framenum]:
                continue

            else:
                gvals,crossings,surr = frame_cont_ent_data[nc]
                Ncrossings, Ccrossings = crossings
                Nsurr, Csurr = surr
                #print(nc, frame_cont_ent_data[nc])
                #print(Ncrossings, Ccrossings)
                #print(Nsurr, Csurr)

                chng_data = []
                #print(nc)
                if nc in frame_Nterm_change_results:
                    chng_data.append(frame_Nterm_change_results[nc])
                    #get frame crossing if chng = 0,1,4
                    #else get ref crossing if chnge = 2,3,
                    #print(ref_cont_ent_data[nc])
                    #print(frame_cont_ent_data[nc])
                    if frame_Nterm_change_results[nc][0] in [0,1,4]:
                        Ncrossings = frame_cont_ent_data[nc][1][0]
                        Nsurr = frame_cont_ent_data[nc][2][0]
                        #print('N-GAIN',nc, frame_Nterm_change_results[nc], Ncrossings)

                    elif frame_Nterm_change_results[nc][0] in [2,3]:
                        Ncrossings = ref_cont_ent_data[nc][1][0]
                        Nsurr = ref_cont_ent_data[nc][2][0]
                        #print('N-LOSS',nc, frame_Nterm_change_results[nc], Ncrossings)

                else:
                    Ncrossings = []
                    Nsurr = []
                    chng_data.append([])

                if nc in frame_Cterm_change_results:
                    chng_data.append(frame_Cterm_change_results[nc])

                    #get frame crossing if chng = 0,1,4
                    #else get ref crossing if chnge = 2,3,
                    if frame_Cterm_change_results[nc][0] in [0,1,4]:
                        Ccrossings = frame_cont_ent_data[nc][1][1]
                        Csurr = frame_cont_ent_data[nc][2][1]
                        #print('C-GAIN',nc, frame_Cterm_change_results[nc],frame_cont_ent_data[nc] ,Ccrossings)

                    elif frame_Cterm_change_results[nc][0] in [2,3]:
                        Ccrossings = ref_cont_ent_data[nc][1][1]
                        Csurr = ref_cont_ent_data[nc][2][1]
                        #print('C-LOSS',nc, frame_Cterm_change_results[nc], frame_cont_ent_data[nc],Ccrossings)

                else:
                    Ccrossings = []
                    Csurr = []
                    chng_data.append([])


                #print(framenum, nc, gvals, chng_data, crossing_data, ent_data)
                outdata[framenum][nc] = [gvals, chng_data, [Ncrossings, Ccrossings], [Nsurr, Csurr]]
                #print(framenum, nc, outdata[framenum][nc])

                #print(outdata[framenum][nc])

        frame_times.append(time.time() - frame_start_time)

    traj_anal_time=time.time()-start_time
    print(f'traj_anal time: {traj_anal_time}')
    print('\n########################################END analysis of trajectory########################################\n')

    mean_frame_time = np.mean(frame_times)
    print(f'mean_frame_time: {mean_frame_time}')


    GQ_outfilename =  f'{out_path}{script_title}_{script_version}/output/{outfile_basename}_GQ.txt'
    np.savetxt(GQ_outfilename, GQ_list, header='framenum, frametime, G, Q')
    print(f'Saved: {GQ_outfilename}')

#output
outfilename = f'{out_path}{script_title}_{script_version}/output/{outfile_basename}.pkl'
with open(outfilename, "wb") as fh:
    pickle.dump(outdata, fh)

print(f'Saved: {outfilename}')

print('\n########################################END output########################################\n')
######################################################################################################################
comp_time=time.time()-start_time
print(f'computation time: {comp_time}')
print(f'NORMAL TERMINATION')
