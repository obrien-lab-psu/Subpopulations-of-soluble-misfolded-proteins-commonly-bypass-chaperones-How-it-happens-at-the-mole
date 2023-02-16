# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:22:02 2021

@author: nissley
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
from matplotlib.pyplot import figure
figure(figsize=(14, 12), dpi=300)

# the fitting function - note that we constrain the two coefficient to sum to unity
def func(t, k1, k2, f1):
    return (f1*np.exp(-k1*t))+((1.0-f1)*np.exp(-k2*t))

ftemp = open('z2.dat')
temp = ftemp.readlines()
ftemp.close()

fit_dict = {}

j = 1
for line in temp:
    junk = line.split()
    fit = []
    for k in junk:
        fit.append(float(k))
    fit_dict[j] = fit
    j += 1

fsize=14
data_color = 'dodgerblue'
fit_color = 'green'
for i in range (1, 20+1):
    
    # get the data for this subplot
    raw = np.loadtxt(str(i)+'.dat')
    time = raw[:,0]
    prob = (100. - raw[:,1])/100.

    # do the fitting
    """
     do the fitting
    sigma = np.ones(len(time))
    sigma[0] = 0.01 # gives higher weights to first point, so it will go almost exactly through it
    popt, pcov = curve_fit(func, time, prob, p0=[1E-5, 1E-8, 0.2], bounds=(0., [np.inf, np.inf, 1.]))
    k1 = popt[0]
    k2 = popt[1]
    f1 = popt[2]
    f2 = 1.0 - popt[2]
    print ('PLOT', i)
    print ('f1:', f1)
    print ('k1:', k1)
    print ('t1:', 1./k1)
    print ('f2:', f2)
    print ('k2:', k2)
    print ('t2:', 1./k2)
    print ()
    fit = (f1*np.exp(-k1*time))+((f2)*np.exp(-k2*time))
    r = pearsonr(fit, prob)[0]
    print ('fit', fit)
    continue
    """
    # make the plot(s)
    ax = plt.subplot(5, 5, int(i))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if i == 1:
        plt.scatter(time, prob, color=data_color, alpha=0.7, label='Experiment')
        plt.scatter(time, fit_dict[i], color=fit_color, alpha=0.7, marker='^', label='Fit to Eq. 1')
        plt.legend(framealpha=0.5, loc='upper right', fancybox=True)
    else:
        plt.scatter(time, prob, color=data_color, alpha=0.7)
        plt.scatter(time, fit_dict[i], color=fit_color, alpha=0.7, marker='^')
    plt.ylim(-0.05, 1.05)
    
    # xticks
    if i == 1:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 10, 20, 30], fontsize=fsize)
    elif i == 2:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60, 80], fontsize=fsize)
    elif i == 3:
        plt.xlim(-1.0, time[-1]+5)
        plt.xticks([0, 50, 100, 150, 200, 250], fontsize=fsize)
    elif i == 4:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 50, 100, 150, 200], fontsize=fsize)
    elif i == 5:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 50, 100, 150], fontsize=fsize)
    elif i == 6:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 50, 100, 150, 200], fontsize=fsize)
    elif i == 7:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 10, 20, 30, 40], fontsize=fsize)
    elif i == 8:
        plt.xlim(-5.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60], fontsize=fsize)
    elif i == 9:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 10, 20, 30, 40, 50], fontsize=fsize)
    elif i == 10:
        plt.xlim(-5.0, time[-1]+4)
        plt.xticks([0, 10, 20, 30, 40, 50, 60], fontsize=fsize)
    elif i == 11:
        plt.xlim(-1.0, time[-1]+1)
        plt.xticks([0, 2, 4, 6, 8, 10], fontsize=fsize)
    elif i == 12:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 10, 20, 30, 40], fontsize=fsize)
    elif i == 13:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60], fontsize=fsize)
    elif i == 14:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60, 80], fontsize=fsize)
    elif i == 15:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60], fontsize=fsize)
    elif i == 16:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60, 80], fontsize=fsize)
    elif i == 17:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 100, 200, 300, 400], fontsize=fsize)
    elif i == 18:
        plt.xlim(-1.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60], fontsize=fsize)
    elif i == 19:
        plt.xlim(-10.0, time[-1]+4)
        plt.xticks([0, 50, 100, 150, 200, 250], fontsize=fsize)
    elif i == 20:
        plt.xlim(-10.0, time[-1]+4)
        plt.xticks([0, 20, 40, 60, 80, 100], fontsize=fsize)
    else:
        print ('Unknown value of i')
        sys.exit()
        
    # axis labels and yticks
    if i == 11:
        plt.xlabel('Time, h', fontsize=fsize)
    #elif i == 12:
     #   plt.xlabel('Time, h', fontsize=fsize)
    else:
        plt.xlabel('Time, min', fontsize=fsize)
    if i in [1, 6, 11, 16]:
        plt.ylabel(r'$P_{\rm NN}(t)$', fontsize=fsize)
        plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=fsize)
    else:
        plt.yticks([0.0, 0.25, 0.50, 0.75, 1.0], [], fontsize=fsize)
plt.tight_layout()
plt.subplots_adjust(hspace=0.5)
plt.savefig('figureS3_v1.png', dpi=300)
