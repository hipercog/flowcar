#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:01:09 2021

@author: bcowley
"""

import csv, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

root = '/media/bcowley/Transcend/flowcar-backup/cogcarsim_physio_2019/sync_physio2ccs/'
ssn = [1, 5, 6, 7, 8]
for s in ssn:
    pth = root + '01/Nexus_data/session_01-0{}'.format(s)
    bios = [f for f in os.listdir(pth) if ('EOG' in f and 'trial' in f)]
    for b in bios:
        fin = open(os.path.join(pth, b))
        eog = list(csv.reader(fin))
        eog = np.array(eog).astype('float')
        
        ts = eog[:, 0]
        eog = eog[:, 1]
        ts_new = np.linspace(ts.min(), ts.max(), 2000) 
        spl = make_interp_spline(ts, eog, k=7)
        eog_smooth = spl(ts_new)
        
        fig = plt.figure(figsize=(80, 60))
        plt.plot(ts_new, eog_smooth)
#        plt.show()
        
        
fin = open(pth + '/' + fnm[0])
bvp = list(csv.reader(fin))
bvpi = np.array(bvp).astype('float')

ts = bvpi[:, 0]

bvp = bvpi[:, 1]

#define x as 200 equally spaced values between the min and max of original x 
ts_new = np.linspace(ts.min(), ts.max(), 2000) 

#define spline with degree k=7
spl = make_interp_spline(ts, bvp, k=7)
bvp_smooth = spl(ts_new)

#create smooth line chart 
#matplotlib(qt)
plt.plot(ts_new, bvp_smooth)
plt.show()