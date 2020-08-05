#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 21:34:44 2020

@author: bcowley
"""

import os, sys, math, json, csv
import numpy as np

#DEBUG
#inpath = '/media/bcowley/Transcend/flowcar-backup/01/Nexus_data/session_01-01/'
#nxfile = 'nexus_1573639466.0630298.txt'

def getframes(YEAR = '2017'): # change to a parameter when ready for both data types
    frames = open(os.path.join(os.path.dirname(__file__), 'Frames' + YEAR + '.csv'))
    frames = list(csv.reader(frames))
    frames = np.asarray(frames)
    return frames

# Get signal out of JSON file
def getsig(nxfile, sigs):
    nxdata = []
    #    create default data in case there are reading problems
    nx_ts_i = {'ts':math.nan}
    nx_sg_i = dict((k,math.nan) for k in sigs)
    with open(nxfile) as f:
        i = 0
        for line in f:
            i += 1
            try:
                json_object = json.loads(line)
                nx_ts_i = json_object[0]
                nx_sg_i = json_object[1]
            except:
                print("JSON line-" + str(i) + " read error for line: " + line)
            else:
                fields = [nx_ts_i[1]]
                for s in sigs:
                    fields.append(nx_sg_i[s])
                nxdata.append(fields)
            
        f.close()
    return nxdata


# Function parses given nexus file, 
def parse_nxs(inf, chans = ['A', 'B', 'E', 'F']):
    dat = getsig(inf, chans)
# parse inf for subject, session, run indices
    
# Get correct frame for this subject, session, run from frames;
    idx = frames[:,0]
    
# Get world time stamp data from pupil folder
# Get the run times from cogcarsim database
# zero the dat timestamp data, and interp based on cogcarsim
    return dat

# Function walks through folders under given path & processes any nexus files found
def traversi(inpath, sigs = ['EDA', 'BVP', 'EOG']):
    sigmap = {'EDA':'E', 'EOG':['A', 'B'], 'BVP':'F'}
    inpath = os.path.abspath(inpath)
    for root, subdirs, files in os.walk(inpath):
        for f in files:
            if f.startswith('nexus') & f.endswith('.txt'):
                nxf = os.path.join(root, f)
                print('Found Nexus, parsing ' + sigs + ' from ' + nxf)
                dat = parse_nxs(nxf, sigmap(sigs))
                for s in sigs:
                    outf = os.path.join(root, s + '_' + f[0:-4] + '.csv')
                    outf = open(outf, "w", newline="")
                    writer = csv.writer(outf)
#                    Extract + write timestamps and correct columns from dat
                    writer.writerows(dat)

            if subdirs != []:
                print('Going deeper')
                for s in subdirs:
                    traversi(s)

#if __name__ == "__main__":
#    parse_nxs(sys.argv[1], sys.argv[2])

#DEBUG
frames = getframes()
#parse_nxs(inpath, 'EOG')
