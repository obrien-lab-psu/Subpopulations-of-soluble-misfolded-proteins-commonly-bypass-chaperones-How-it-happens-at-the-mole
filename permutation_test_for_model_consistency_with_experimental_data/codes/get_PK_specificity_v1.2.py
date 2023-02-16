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

np.set_printoptions(precision=3)

if len(sys.argv) != 3:
    print('[1] path to stephens data')
    print('[2] outpath')
    sys.exit()

def load_LipMS_data(file_path):
    print('--------------------------------------------------------------------------------------------------------------')
    print('\nLoading LipMS data')
    files = glob.glob(file_path)
    LipMS_data = []
    for f in files:
        excel_df = pd.read_excel(f)
        excel_df = excel_df[['proteinaseKsite']]

        for peptide in excel_df['proteinaseKsite']:
            if type(peptide) == str:
                if '[' in peptide:
                    pass
                else:
                    AA = "".join(re.split("[^a-zA-Z]*",peptide))

                    LipMS_data += [AA]

    return LipMS_data


#load user arguments
global tart, end, aID
outpath = sys.argv[2]

#load lipms data
LipMS_data = load_LipMS_data(sys.argv[1])

AA_u, AA_counts = np.unique(LipMS_data, return_counts = True)
AA_prob = AA_counts/np.sum(AA_counts)

ofile = f'{outpath}'
np.savetxt(ofile, np.hstack((AA_u[:,None], AA_counts[:,None], AA_prob[:,None])), header=f'AA, count, prob', fmt='%s')
print(f'SAVED: {ofile}')

print('NORAML TERMINATION')
