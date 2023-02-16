#!/usr/bin/env python3

import numpy as np
import sys,os
import glob
import pickle
import pandas as pd
import re
from scipy.spatial import distance
from itertools import product
import time
import re
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)
pd.options.display.max_columns = 2000

if len(sys.argv) != 3:
    print('[1] path to stephens data')
    print('[2] outfile basename')
    sys.exit()


def load_sig_LipMS_data(file_path):
    print('--------------------------------------------------------------------------------------------------------------')
    print('\nLoading sig LipMS data')
    files = glob.glob(file_path)
    LipMS_data = {}
    for f in files:
        print(f)
        excel_df = pd.read_excel(f)
        LipMS_data[f] = []

        for i in excel_df.index:
            PKsite =  excel_df['proteinaseKsite'][i]
            if '[' in PKsite:
                continue

            peptide = excel_df['Peptide Sequence'][i]
            peptide = "".join(re.findall("[a-zA-Z]+", peptide))

            num_RK = len([x for x in peptide if x == 'R' or x == 'K'])
            print(peptide, PKsite, len(peptide), num_RK)

            LipMS_data[f] += [[len(peptide), num_RK]]


    return LipMS_data

def edges(start, width, end):
    out = []
    i = start
    while i <= end:
        out += [i]
        i += width

    return np.asarray(out)

############
### MAIN ###
############
start_time = time.time()

#load user arguments
outfile = sys.argv[2]

#load lipms data
LipMS_data = load_sig_LipMS_data(sys.argv[1])
print(LipMS_data)
outdata = []

for k,v in LipMS_data.items():
    print(k)
    print(v)
    outdata.append(v)

out = {}
outdata = np.vstack(outdata)
print(min(outdata[:,0]), max(outdata[:,0]))
print(min(outdata[:,1]), max(outdata[:,1]))

xedges = edges(0.5, 1, 50.5)
print(xedges)
yedges = edges(0.5, 1, 5.5)
print(yedges)
H, xedges, yedges = np.histogram2d(outdata[:,0], outdata[:,1], bins=[xedges, yedges])
print(H)
P = H/np.sum(H)
print(P)
out['H'] = H
out['P'] = P
out['xedges'] = xedges
out['yedges'] = yedges

outfilename = f'{outfile}.pkl'
with open(outfilename, "wb") as fh:
        pickle.dump(out, fh)

print(f'Saved: {outfilename}')
end_time = time.time() - start_time
print(f'NORAML TERMINATION: {end_time}')
