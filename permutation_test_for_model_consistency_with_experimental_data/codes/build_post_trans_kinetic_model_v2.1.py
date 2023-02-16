#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import pyemma as pem
import parmed as pmd
import mdtraj as mdt
import msmtools
import time

################################# Arguments ###################################
start_time = time.time()
# Default values
end_t = 60 # in seconds
dt = 0.015/1000
nsave = 5000
alpha = 4331293.0
n_window = 200
n_traj = 100
ens_label_list = ['wt']
n_cluster = 400
stride=10
n_large_states = 10
n_small_states = 2
lag_t = 1
start_idx = 1
end_idx = 1
sample_size = 5
dcd_dir = None
psf_file = None
cor_file = None
GQ_ent_dir = None
if_cluster = True
if_sample = True
start_frame = 0
end_frame = 9999

# read control file
ctrlfile = ''

if len(sys.argv) == 1:
    print('[1] = control file')
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:", ["ctrlfile="])
except getopt.GetoptError:
    print('[1] = control file')
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print('[1] = control file')
        sys.exit()
    elif opt in ("-f", "--ctrlfile"):
        ctrlfile = arg

if not os.path.exists(ctrlfile):
    print('Error: cannot find control file ' + ctrlfile + '.')
    sys.exit()

file_object = open(ctrlfile,'r')
try:
    for line in file_object:
        line = line.strip()
        if not line:
            # This is a blank line
            continue
        if line.startswith('#'):
            # This is a comment line
            continue
        if line.startswith('end_t'):
            words = line.split('=')
            end_t = float(words[1].strip())
            continue
        if line.startswith('dt'):
            words = line.split('=')
            dt = float(words[1].strip())
            continue
        if line.startswith('nsave'):
            words = line.split('=')
            nsave = int(words[1].strip())
            continue
        if line.startswith('alpha'):
            words = line.split('=')
            alpha = float(words[1].strip())
            continue
        if line.startswith('n_window'):
            words = line.split('=')
            n_window = int(words[1].strip())
            continue
        if line.startswith('n_traj'):
            words = line.split('=')
            n_traj = int(words[1].strip())
            continue
        if line.startswith('ens_label_list'):
            words = line.split('=')
            ens_label_list = words[1].strip().split()
            continue
        if line.startswith('n_cluster'):
            words = line.split('=')
            n_cluster = int(words[1].strip())
            continue
        if line.startswith('stride'):
            words = line.split('=')
            stride = int(words[1].strip())
            continue
        if line.startswith('n_large_states'):
            words = line.split('=')
            n_large_states = int(words[1].strip())
            continue
        if line.startswith('n_small_states'):
            words = line.split('=')
            n_small_states = int(words[1].strip())
            continue
        if line.startswith('lag_t'):
            words = line.split('=')
            lag_t = int(words[1].strip())
            continue
        if line.startswith('start_idx'):
            words = line.split('=')
            start_idx = int(words[1].strip())
            continue
        if line.startswith('end_idx'):
            words = line.split('=')
            end_idx = int(words[1].strip())
            continue
        if line.startswith('sample_size'):
            words = line.split('=')
            sample_size = int(words[1].strip())
            continue
        if line.startswith('start_frame'):
            words = line.split('=')
            start_frame = int(words[1].strip())
            continue
        if line.startswith('end_frame'):
            words = line.split('=')
            end_frame = int(words[1].strip())
            continue
        if line.startswith('dcd_dir'):
            words = line.split('=')
            dcd_dir = words[1].strip().split()
            continue
        if line.startswith('psf_file'):
            words = line.split('=')
            psf_file = words[1].strip()
            continue
        if line.startswith('cor_file'):
            words = line.split('=')
            cor_file = words[1].strip()
            continue
        if line.startswith('GQ_ent_dir'):
            words = line.split('=')
            GQ_ent_dir = words[1].strip()
            continue
        if line.startswith('outpath'):
            words = line.split('=')
            outpath = words[1].strip()
            continue
        if line.startswith('outfile_basename'):
            words = line.split('=')
            outfile_basename = words[1].strip()
            continue
        if line.startswith('if_cluster'):
            words = line.split('=')
            if_cluster = int(words[1].strip())
            if if_cluster == 1:
                if_cluster = True
            elif if_cluster == 0:
                if_cluster = False
            else:
                print('Error: if_cluster can only be either 0 or 1.')
                sys.exit()
            continue
        if line.startswith('if_sample'):
            words = line.split('=')
            if_sample = int(words[1].strip())
            if if_sample == 1:
                if_sample = True
            elif if_sample == 0:
                if_sample = False
            else:
                print('Error: if_sample can only be either 0 or 1.')
                sys.exit()
            continue
finally:
     file_object.close()

dt = dt*nsave*alpha/1e9 # in seconds
################################# Functions ###################################
def standardize(data):
    data_con = data[0]
    for i in range(1, len(data)):
        data_con = np.vstack((data_con, data[i]))
    data_mean = np.mean(data_con, axis=0)
    data_std = np.std(data_con, axis=0)
    result = [(d - data_mean) / data_std for d in data]
    return [result, data_mean, data_std]

def unstandardize(data, data_mean, data_std):
    result = data * data_std + data_mean
    return result

# Building Master Equations
def dP(t, P, M):
    return np.dot(M,P)

def Master_equation(t_span, K, P0, t_eval):
    sol = solve_ivp(dP, t_span, P0, t_eval=t_eval, args=(K,))
    return sol

def estimate_rate_matrix(ens_label_type, dtrajs, n_state, dt):
    ksum_matrix = np.zeros((n_state, n_state))
    P_matrix = np.zeros((n_state, n_state))
    C_matrix = msmtools.estimation.count_matrix(dtrajs, 1)
    C_matrix = C_matrix.toarray()
    if len(C_matrix) != n_state:
        for i in range(len(C_matrix), n_state):
            C_matrix.append(np.zeros((1,C_matrix.shape[1])), axis=0)
            C_matrix.append(np.zeros((C_matrix.shape[0],1)), axis=1)
    for i in range(n_state):
        if C_matrix[i,i] != 0:
            ksum_matrix[i,i] = -math.log(C_matrix[i,i]/np.sum(C_matrix[i,:]))/dt
        else:
            ksum_matrix[i,i] = 2/dt # assume the mean dwell time of this state is dt/2
        for j in range(n_state):
            if i != j:
                P_matrix[j,i] = C_matrix[i,j]
        if np.sum(P_matrix[:,i]) != 0:
            P_matrix[:,i] /= np.sum(P_matrix[:,i])

    rate_matrix = np.dot(P_matrix, ksum_matrix)
    # enforce native state to be a sink
    # for i in range(n_state):
    #     rate_matrix[i,-1] = 0
    for i in range(n_state):
        rate_matrix[i,i] = -np.sum(rate_matrix[:,i])

    return rate_matrix

def exp_fun(x, k):
    return np.exp(-k*x)


################################## MAIN #######################################
if not if_cluster:
    npzfile = np.load('msm_data.npz', allow_pickle=True)

prefix_file = '_act'

G_list_0_list = []
cor_list = []
trajid_list = []
mtype2trajid = []
rate_matrix_list = []
xlim_list = []
ylim_list = []
fig_list = []
fig_tpt_list = []
traj_idx_2_mtype_idx = {}

# combine trajs and do clustering and PCCA++
for i_ax, ens_label_type in enumerate(ens_label_list):
    print(i_ax, ens_label_type)
    max_length = len(pmd.load_file(psf_file).atoms)

    files = os.listdir(GQ_ent_dir)
    files = sorted([x for x in os.listdir(GQ_ent_dir) if x.endswith('_GQ.txt') and ens_label_type in x], key=lambda x: int(x.split('_')[3].strip('t')))
    files = [f'{GQ_ent_dir}{f}' for f in files]

    for i in range(len(files)):
        i += len(cor_list)
        traj_idx_2_mtype_idx[i] = i_ax

    mtype2trajid.append([i+len(cor_list) for i in range(len(files))])

    for fi, f in enumerate(files):
        print(f)
        cor_list += [np.loadtxt(f)[start_frame:end_frame,[3,2]]]
    trajid_list += [i for i in range(len(files))]

print(f'trajid_list: {trajid_list}')
print(f'mtype2trajid: {mtype2trajid}')
print(len(cor_list))

#qualtiy check cor_list for lengths of data
print(f'start_frame: {start_frame} | end_frame: {end_frame}')
for data in cor_list:
    print(data.shape)

#Clustering
if if_cluster:
    std_cor_list, cor_mean, cor_std = standardize(cor_list)
    cluster = pem.coordinates.cluster_kmeans(std_cor_list, k=n_cluster, max_iter=5000, stride=stride)
    dtrajs = cluster.dtrajs
    center = unstandardize(cluster.clustercenters, cor_mean, cor_std)
else:
    dtrajs = list(npzfile['dtrajs'])
    center = npzfile['center']

# Get connective groups and build MSMs
c_matrix = msmtools.estimation.count_matrix(dtrajs, lag_t).toarray()
sub_groups = msmtools.estimation.connected_sets(c_matrix)
active_groups = []
for sg in sub_groups:
    for ssg in sg:
        tag_found = False
        for dtraj in dtrajs:
            if ssg in dtraj:
                tag_found = True
                break
        if not tag_found:
            break
    if tag_found:
        active_groups.append(sg)
print('Total number of active groups: %d'%(len(active_groups)))

msm_list = []
for ag in active_groups:
    cm = msmtools.estimation.largest_connected_submatrix(c_matrix, lcc=ag)
    if len(cm) == 1:
        msm = None
    else:
        T = msmtools.estimation.transition_matrix(cm, reversible=True)
        msm = pem.msm.markov_model(T, dt_model=str(dt)+' s')
    msm_list.append(msm)

meta_dist = []
meta_set = []
eigenvalues_list = []
for idx_msm, msm in enumerate(msm_list):
    if idx_msm == 0:
        n_states = n_large_states
    else:
        n_states = n_small_states
    if msm == None:
        eigenvalues_list.append(None)
        dist = np.zeros(n_cluster)
        iidx = active_groups[idx_msm][0]
        dist[iidx] = 1.0
        meta_dist.append(dist)
        meta_set.append(active_groups[idx_msm])
    else:
        eigenvalues_list.append(msm.eigenvalues())
        # coarse-graining
        while n_states > 1:
            tag_empty = False
            pcca = msm.pcca(n_states)
            for ms in msm.metastable_sets:
                if ms.size == 0:
                    tag_empty = True
                    break
            if not tag_empty:
                break
            else:
                n_states -= 1
                print('Reduced number of states to %d for active group %d'%(n_states, idx_msm+1))
        if n_states == 1:
            # use observation prob distribution for non-active set
            dist = np.zeros(n_cluster)
            for nas in active_groups[idx_msm]:
                for dtraj in dtrajs:
                    dist[nas] += np.count_nonzero(dtraj == nas)
            dist /= np.sum(dist)
            meta_dist.append(dist)
            meta_set.append(active_groups[idx_msm])
        else:
            for i, md in enumerate(msm.metastable_distributions):
                dist = np.zeros(n_cluster)
                s = np.sum(md[msm.metastable_sets[i]])
                set_0 = []
                for idx in msm.metastable_sets[i]:
                    iidx = active_groups[idx_msm][idx]
                    dist[iidx] = md[idx]
                    set_0.append(iidx)
                dist = dist / s
                meta_dist.append(dist)
                meta_set.append(set_0)
meta_dist = np.array(meta_dist)
meta_set = np.array(meta_set)

coarse_state_centers = center[meta_dist.argmax(1)]
cg_center_order_idx = np.argsort(coarse_state_centers[:,0])
micro_to_meta = np.zeros(n_cluster)
meta_set = meta_set[cg_center_order_idx]
meta_dist = meta_dist[cg_center_order_idx, :]
for idx, ms in enumerate(meta_set):
    for mms in ms:
        micro_to_meta[mms] = idx
meta_dtrajs = []
for traj in dtrajs:
    meta_traj = np.zeros(len(traj), dtype=int)
    for i, idx in enumerate(traj):
        meta_traj[i] = micro_to_meta[idx]
    meta_dtrajs.append(meta_traj)
meta_dtrajs = np.array(meta_dtrajs)

n_states = len(meta_set)
print('Total %d metastable states were grouped'%n_states)

if if_sample:
    cluster_indexes = pem.util.discrete_trajectories.index_states(dtrajs)
    if len(cluster_indexes) < n_cluster:
        cluster_indexes = list(cluster_indexes)
        for i in range(len(cluster_indexes), n_cluster):
            cluster_indexes.append(np.array([[]]))
        cluster_indexes = np.array(cluster_indexes)
    samples = pem.util.discrete_trajectories.sample_indexes_by_distribution(cluster_indexes,
                                                                            meta_dist,
                                                                            sample_size)
else:
    samples = npzfile['meta_samples']
meta_samples = samples

sampled_traj = None
print(dcd_dir)
#dcd_files = os.listdir(f'{dcd_dir}')
#print(dcd_files)
for i, meta_state in enumerate(samples):
    print('\nGet samples of metastable state %d'%(i+1))
    print('msm_traj_idx, traj_idx, frame_idx, loc_frame_idx, traj_idx_1, label')
    for idx in meta_state:
        msm_traj_idx = idx[0]
        frame_idx = idx[1]

        loc_frame_idx = start_frame + frame_idx
        traj_idx = trajid_list[msm_traj_idx]
        traj_idx_1 = traj_idx_2_mtype_idx[msm_traj_idx]
        #traj_idx_1 = int(traj_idx / (end_idx-start_idx+1))
        #traj_idx_2 = int(traj_idx - traj_idx_1 * (end_idx-start_idx+1))
        label = ens_label_list[traj_idx_1]
        print(msm_traj_idx, traj_idx, frame_idx, loc_frame_idx, traj_idx_1,  label)
        dcd_files = os.listdir(f'{dcd_dir[traj_idx_1]}')
        traj_file_path = [f'{dcd_dir[traj_idx_1]}{dcd}' for dcd in dcd_files if dcd.startswith(f'{traj_idx+1}_') and label in dcd and dcd.endswith('.dcd')][0]
        print(traj_file_path)

        if traj_file_path == '':
            print(f'DCD not found for {label} traj# {traj_idx_1}')
            quit()
        traj = mdt.load(traj_file_path, top=psf_file)
        sel = traj.topology.select('resid 0 to '+str(max_length-1))
        traj = traj.atom_slice(sel)
        if sampled_traj is None:
            sampled_traj = traj[loc_frame_idx]
        else:
            sampled_traj += traj[loc_frame_idx]

sampled_traj = sampled_traj.center_coordinates()
sampled_traj = sampled_traj.superpose(sampled_traj)
sampled_traj.save(f'{outpath}{outfile_basename}_sampled_traj.dcd', force_overwrite=True)
print(f'\nSAVED: {outpath}{outfile_basename}_sampled_traj.dcd')

if_entangled_list = [False for i in range(n_states)]
# analysis MSM for each ens_label
for i_ax, ens_label_type in enumerate(ens_label_list):

    state_indices = []
    for state_id in range(0, n_states):
        state_indices.append([])
        for i_1, md in enumerate(meta_dtrajs[mtype2trajid[i_ax]]):
            for i_2, mdd in enumerate(md):
                if mdd == state_id:
                    state_indices[-1].append([mtype2trajid[i_ax][i_1], i_2])

    max_T_len = 0
    for i in mtype2trajid[i_ax]:
        if meta_dtrajs[i].shape[0] > max_T_len:
            max_T_len = meta_dtrajs[i].shape[0]

    PPT = np.zeros((max_T_len, n_states))
    meta_dtrajs_extended = []
    for i in mtype2trajid[i_ax]:
        (N, be) = np.histogram(meta_dtrajs[i][-n_window:], bins=np.arange(-0.5, n_states, 1))
        meta_dtraj_last = np.argwhere(N == np.max(N))[0][0]
        for j in range(max_T_len):
            if j >= len(meta_dtrajs[i]):
                state_0 = meta_dtraj_last
            else:
                state_0 = meta_dtrajs[i][j]
            PPT[j,state_0] += 1
        mde = []
        for j in range(max_T_len):
            if j >= len(meta_dtrajs[i]):
                state_0 = meta_dtraj_last
            else:
                state_0 = meta_dtrajs[i][j]
            mde.append(state_0)
        meta_dtrajs_extended.append(mde)

    rate_matrix = estimate_rate_matrix(ens_label_type, meta_dtrajs_extended, n_states, dt)
    rate_matrix_list.append(rate_matrix)

    fff = open(f'{outpath}{ens_label_type}_MSTS.dat', 'w')
    for i in range(PPT.shape[0]):
        PPT[i,:] = PPT[i,:] / np.sum(PPT[i,:])
        fff.write('%8.4f '%((i+1)*dt))
        for j in range(n_states):
            fff.write('%6.4f '%PPT[i,j])
        fff.write('\n')
    fff.close()
    print(f'SAVED: {outpath}{ens_label_type}_MSTS.dat')

    # build master equation
    t_span = [0, end_t]
    t_eval = np.linspace(0, end_t, 1000)
    P0 = PPT[0,:]
    sol = Master_equation(t_span, rate_matrix, P0, t_eval)
    t_span = sol.t
    PP = sol.y

    fff = open(f'{outpath}{ens_label_type}_METS.dat', 'w')
    for i in range(t_span.shape[0]):
        fff.write('%8.4f '%(t_span[i]))
        for j in range(n_states):
            fff.write('%6.4f '%abs(PP[j,i]))
        fff.write('\n')
    fff.close()
    print(f'SAVED: {outpath}{ens_label_type}_METS.dat')

np.savez(f'{outpath}{outfile_basename}_msm_data.npz',
         dtrajs = dtrajs,
         center = center,
         eigenvalues_list = eigenvalues_list,
         meta_dtrajs = meta_dtrajs,
         meta_set = meta_set,
         meta_samples = meta_samples,
         rate_matrix_list = rate_matrix_list)
print(f'SAVED: {outpath}{outfile_basename}_msm_data.npz')
print(f'\nNORMAL TERMINATION @ {time.time() - start_time}')
