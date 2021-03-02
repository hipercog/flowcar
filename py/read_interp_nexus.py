#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 21:34:44 2020

@author: bcowley
"""

import os, math, json, csv, sys
import numpy as np
import logging as lg
import datetime as dt
import matplotlib.pyplot as plt
import scipy.interpolate as spi



# Load the manually-observed frames information
def getframes(inp, YEAR = '2019'): 
    frames = open(os.path.join(inp, 'frames' + YEAR + '.csv'))
    frames = list(csv.reader(frames))
    frames = np.asarray(frames)
    
    BLTs = open(os.path.join(inp, 'baselineTimes' + YEAR + '.csv'))
    BLTs = list(csv.reader(BLTs))
    BLTs = np.asarray(BLTs)
    
    return frames, BLTs


# Load the database tables for CogCarSim logs 
#  - these must have already been exported from SQL to csv
def getCogCarSimData(inp):
    CCSsteps = np.array(np.genfromtxt(os.path.join(inp, 'Databases'
              , 'cogcarsim2_2019_step.csv'), delimiter=",", skip_header = 1))
    CCSruns = open(os.path.join(inp, 'Databases', 'cogcarsim2_2019_run.csv'))
    CCSruns = list(csv.reader(CCSruns))
    CCSruns = np.asarray(CCSruns)
    return CCSruns, CCSsteps



# DECLARE GLOBALS
inpath = '/media/bcowley/Transcend/flowcar-backup/'
outpath = os.path.join(inpath, 'read_interp_nexus')
sigmap = {'EDA':'E', 'EOG':['A', 'B'], 'BVP':'F'}
mapsig = {'E':'EDA', 'A':'EOG', 'B':'EOG', 'F':'BVP'}
siglist = ['EDA', 'BVP', 'EOG']
if "frames" not in globals() and "BLTs" not in globals():
    frames, BLTs = getframes(inpath)
if "CCSrun" not in globals() and "CCSstp" not in globals():
    CCSrun, CCSstp = getCogCarSimData(inpath)



# Get signal out of JSON file
def getsig(nexusfile, chans = ['A', 'B', 'E', 'F']):
    nxdata = []
    #    create default data in case there are reading problems
    nx_ts_i = {'ts':math.nan}
    nx_sg_i = dict((k,math.nan) for k in chans)
    with open(nexusfile) as f:
        i = 0
        for line in f:
            i += 1
            try:
                json_object = json.loads(line)
                nx_ts_i = json_object[0]
                nx_sg_i = json_object[1]
            except:
                lg.exception("JSON read error line{}: {}".format(str(i), line))
            else:
                fields = [nx_ts_i[1]]
                for c in chans:
                    fields.append(nx_sg_i[c])
                nxdata.append(fields)
            
        f.close()
        
    nxarr = np.array(nxdata).T
    idx = np.divide([len(np.unique(r)) * 100 for r in nxarr], len(nxarr.T)) > 5
        
    return np.array(chans)[idx[1:]], nxarr[idx,:]


# Get Pupil session info
def getpupilsystime(infoplyrjsonfile):
    
    info = json.load(open(infoplyrjsonfile))
    
    for k, v in info.items():
        if k == "start_time_synced_s":
            sync = v
        elif k == "start_time_system_s":
            syst = v
    
    systime = syst - sync
    
    return systime


# Function to plot and save figure comparing pre- and post-interpolation data
def sanity_check(svpath, signame, t1, d1, t2, d2):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    
    ax1.plot(t1, 1.0/d1, '-')
    ax1.set_title('raw ' + signame)
    
    ax2.plot(t2, 1.0/d2, '-')
    ax2.get_yaxis().set_visible(False)
    ax2.set_title('interpolated ' + signame)
    
    if not os.path.exists(os.path.dirname(svpath)):
        os.makedirs(os.path.dirname(svpath), exist_ok = True)
        
    plt.savefig(svpath + '_interp_' + signame)
    plt.show()


# Function parses given nexus file, 
def parse_nxs(inp, fnm, sigs):
    
    # parse fnm for subject, session, run indices
    pfx = os.sep + 'Nexus_data' + os.sep + 'session_'
    pfx = inp.find(pfx) + len(pfx)
    sbj = inp[pfx:pfx + 2]
    ssn = inp[pfx + 3:pfx + 5]
    
    lg.info('Subject: ' + sbj)
    lg.info('Session: ' + ssn)
    
    # Get world time stamp data from pupil folder
    pupilpath = os.path.join(inp.replace('Nexus', 'Pupil'), '000')
    world_times = np.load(os.path.join(pupilpath, "world_timestamps.npy"))
    world_times += getpupilsystime(os.path.join(pupilpath, "info.player.json"))

    # get the signal data    
    chx = []
    for k, v in sigmap.items():
        if k in sigs:
            [chx.append(c) for c in v]

    hdr, dat = getsig(os.path.join(inp, fnm), chx)
    
    # Get baseline segment based on BLTs
    frblt = [int(BLTs[np.where(BLTs[:,0] == sbj + '-' + ssn[1] + '-S'), 1]), 
             int(BLTs[np.where(BLTs[:,0] == sbj + '-' + ssn[1] + '-E'), 1])]
    if len(frblt) == 0:
        lg.error('No baseline frames found for subject {}, session {}'.format(
                sbj, ssn))
    bldat = dat[:, frblt[0]:frblt[1]]
    bldat[0, :] -= bldat[0, 0] # zero the timestamps

    # Split time and data to make interpolation of the trials readable
    ts = dat[0,:]
    dat = np.delete(dat, 0, 0)

    # Get correct frame/timebase for this subject, session, & run
    frssn = [len(world_times) - 1]
    tbase = [world_times[-1]]
    for f in np.flip(frames[:,0]):
        if (f.startswith(sbj + '-' + ssn[1] + '-')):
            frssn = [int(frames[np.where(frames[:,0] == f),1])] + frssn
            tbase = [world_times[frssn[0]]] + tbase
    if len(frssn) == 1:
        lg.error('No trial frames found for subject {}, session {}'.format(
                sbj, ssn))
    
    outdat = [None] * (len(frssn) - 1)
    
    # Step through each run and interpolate each signal, building output data
    for i in range(len(frssn) - 1):
        lg.info('Trial {}: frame {} to frame {}'.format(
                i + 1, frssn[i], frssn[i + 1]))
        
        # Get the run times from cogcarsim database
        runidx = CCSrun[:,1] == sbj + '-0' + ssn[1] + '-0' + str(i + 1)
        runTS = CCSstp[CCSstp[:, 1] == int(CCSrun[runidx,0]), 2].ravel()
        runTS -= runTS[0]
        lg.info("Length of run: {} seconds".format(round(runTS[-1])))
        
        # Extract biosigs and physio timestamps for trial i
        t1 = np.argmin(abs(ts - tbase[i]))
        t2 = np.argmin(abs(ts - tbase[i + 1]))
        fi_ts = ts[t1 - 128:t2] - tbase[i]
        fi_dat = dat[:, t1 - 128:t2]
        
        framei = [None] * (len(hdr) + 1)
        framei[0] = runTS

        # extract the right physio data, and interp based on cogcarsim
        for chx in range(len(hdr)):
            lg.info("Extracting {}".format(mapsig[hdr[chx]]))
            #interpolate biosig & get final signal with cogcarsim timestamps
            interper = spi.interp1d(fi_ts, fi_dat[chx,:], kind = 1)
            framei[chx + 1] = interper(runTS)
            # plot comparison of pre- to post-interp
            sanity_check(os.path.join(inp, 'figs', 'trial_{}'.format(i + 1))
              , mapsig[hdr[chx]], fi_ts, fi_dat[chx,:], runTS, framei[chx + 1])
            
        outdat[i] = np.array(framei)
    
    return hdr, outdat, bldat


# Function walks through folders under given path & processes any nexus files
def traversi(indir, outdir, sigs = ['EDA', 'BVP', 'EOG']):
    indir = os.path.abspath(indir)
    for root, subdirs, files in os.walk(indir):
        for f in files:
            if f.startswith('nexus') & f.endswith('.txt'):
                lg.info('Read {} from {}'.format(sigs, os.path.join(root, f)))
                
                # Manage file writing paths etc
                if not os.path.exists(root.replace(indir, outdir)):
                    os.makedirs(root.replace(indir, outdir), exist_ok = True)
                lg.info('Write to {}'.format(root.replace(indir, outdir)))
                ts = dt.datetime.fromtimestamp(
                        float(f[f.find('_') + 1:f.find('.txt')]))
                
                # Extract + write timestamps and correct columns
                hdr, dat, bl = parse_nxs(root, f, sigs)
                for ch in hdr:
                    nm = 'Nexus{}_{}_baseline.csv'.format(
                            mapsig[ch], ts.strftime('%Y%m%d'))
                    outf = os.path.join(root.replace(indir, outdir), nm)
                    outf = open(outf, "w", newline="")
                    writer = csv.writer(outf)
                    fi_ch_bl = bl[[True] + list(hdr==ch), :]
                    writer.writerows(fi_ch_bl.T)
                    
                    
                for i in range(len(dat)):
                    fi_dat = dat[i]
                    for ch in hdr:
                        nm = 'Nexus{}_{}_trial-{}.csv'.format(
                                mapsig[ch], ts.strftime('%Y%m%d'), str(i + 1))
                        outf = os.path.join(root.replace(indir, outdir), nm)
                        outf = open(outf, "w", newline="")
                        writer = csv.writer(outf)
                        fi_ch_dat = fi_dat[[True] + list(hdr==ch), :]
                        writer.writerows(fi_ch_dat.T)

        if subdirs != []:
            lg.info('Going deeper')
            for s in subdirs:
                nu_ind = os.path.join(indir, s)
                traversi(nu_ind, nu_ind.replace(indir, outdir))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        inpath = sys.argv[1]
        outpath = os.path.join(inpath, sys.argv[0])
        if len(sys.argv) > 2:
            siglist = sys.argv[2]


lg.basicConfig(filename=os.path.join(inpath, 'read_interp_nexus.log'), filemode
       = 'w', format='%(name)s - %(levelname)s - %(message)s', level=lg.DEBUG)

#DEBUG
DEBUG = True
if DEBUG:
    traversi(os.path.join(inpath, '01/Nexus_data/session_01-01/')
            , os.path.join(inpath, '01/Nexus_data/session_01-01/read_interp_nexus')
            , siglist)
#    parse_nxs(os.path.join(inpath, '01/Nexus_data/session_01-01/')
#            , 'nexus_1573639466.0630298.txt'
#            , siglist)
else:
    traversi(inpath, outpath, siglist)