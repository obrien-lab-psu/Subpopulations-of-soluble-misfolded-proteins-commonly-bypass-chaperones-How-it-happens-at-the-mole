#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback
import numpy as np
from scipy.integrate import solve_ivp
import scipy.stats as stats

def ODE_fun(t, P, K):
    # P is the vector of concentrations
    # P[0] is [G], single-ring GroEL
    # P[1] is [G|7ATP]
    # P[2] is [G|7ATP|U]
    # P[3] is [G|7ATP|U|ES]
    # P[4] is [G|7ADP|U|ES]
    # P[5] is [ES], GroES
    # P[6] is [7ATP]
    # P[7] is [F], folded protein
    # P[8] is [M], misfolded protein
    # P[9] is [U], unfolded protein
    # P[10] is [G|7ADP]
    
    # K is the vector of rates and prefactors
    # K[0] is ATP binding rate
    # K[1] is unfolded protein binding rate
    # K[2] is GroES binding rate
    # K[3] is ATP hydrolysis rate
    # K[4] is protein release rate
    # K[5] is folded propotion in GroEL
    # K[6] is misfolded propotion in GroEL
    # unfolded propotion in GroEL = 1 - K[5] - K[6]
    # K[7] is folding rate in Bulk
    # K[8] is folded propotion in Bulk
    # misfolded propotion in bulk = 1 - K[8]
    # K[9] is basal ATP hydrolysis rate
    # K[10] is ADP release rate 
    
    dP = np.zeros(len(P))
    dP[0] = -K[0]*P[0]*P[6] + K[4]*P[4] + K[10]*P[10]
    dP[1] = K[0]*P[0]*P[6] - K[1]*P[1]*P[9] - K[9]*P[1]
    dP[2] = K[1]*P[1]*P[9] - K[2]*P[2]*P[5]
    dP[3] = K[2]*P[2]*P[5] - K[3]*P[3]
    dP[4] = K[3]*P[3] - K[4]*P[4]
    dP[5] = -K[2]*P[2]*P[5] + K[4]*P[4]
    dP[6] = -K[0]*P[0]*P[6]
    dP[7] = K[5]*K[4]*P[4] + K[8]*K[7]*P[9]
    dP[8] = K[6]*K[4]*P[4] + (1-K[8])*K[7]*P[9]
    dP[9] = -K[1]*P[1]*P[9] + (1-K[5]-K[6])*K[4]*P[4] - K[7]*P[9]
    dP[10] = K[9]*P[1] - K[10]*P[10]
    
    return dP

def Kinetic_model(t_span, K, P0, t_eval):
    sol = solve_ivp(ODE_fun, t_span, P0, t_eval=t_eval, args=(K,), method='LSODA')
    return sol

############ MAIN ################
# rates are in unit of 1/uM/s or 1/s
rates = np.array([0.2, 2, 20, 0.12, 1, 0, 0, 0, 0, 5/60/28, 0.3])*60 # convert to 1/min

# Bulk folding rate in 1/min
rates[7] = 0.15 # Ritaban, change this for a given protein;
# Folded propotion in Bulk
rates[8] = 0.03 # Ritaban, change this for a given protein;

# Experimental time course of non-native protein
data = np.array([[1.47429621, 0.933879137],
                 [4.98946953, 0.728441715],
                 [10.3912668,  0.522870898],
                 [14.9298449,  0.381040335],
                 [30.1770935,  0.252602994],
                 [44.9022560,  0.209108293],
                 [60.1524734,  0.208029270],
                 [90.2049789,  0.224770784]])


t_eval = data[:,0] # in min

# initial concentrations in the unit of uM
P0 = np.zeros(11)
P0[0] = 0.69 # [G]
P0[5] = 0.69 # [ES]
P0[6] = 1000 # [ATP]
P0[9] = 0.46 # [U]

# Single-ring concentration
P0[0] = P0[0]*2
# 7 ATPs bind at the same time (???)
P0[6] = P0[6]/7

f = open('LOG', 'w')
f.close()
# Solve the equations
for phi_F in np.arange(0.001, 1.001-data[-1,-1], 0.001):
#for phi_F in np.arange(0.05, 1.05, 0.05):
    phi_M = data[-1,-1]/(1-data[-1,-1])*phi_F
    #phi_M = 0
    if phi_M + phi_F >= 1:
        break
    K = np.copy(rates)
    K[5] = phi_F
    K[6] = phi_M
    t_span = [0, t_eval[-1]]
    sol = Kinetic_model(t_span, K, P0, t_eval)
    p_uf = 1-sol.y[7,:]/P0[9]
    (r, p) = stats.pearsonr(data[:,1], p_uf)
    err = np.sum(np.abs(data[:,1] - p_uf))
    f = open('LOG', 'a')
    print('Unfolded fraction: ', p_uf)
    print('[ATP]: ', sol.y[6,:]*7)
    print('phi_F = %f, phi_M = %f, r = %f, error = %f\n'%(phi_F, phi_M, r, err))
    f.write('Unfolded fracion: '+str(p_uf)+'\n')
    f.write('[ATP]: '+str(sol.y[6,:]*7)+'\n')
    f.write('phi_F = %f, phi_M = %f, r = %f, error = %f\n\n'%(phi_F, phi_M, r, err))
    f.close()
        



