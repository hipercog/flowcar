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
    frames = open(os.path.join(inp, 'trialFrames' + YEAR + '.csv'))
    frames = list(csv.reader(frames))
    frames = np.asarray(frames)
    
    BLTs = open(os.path.join(inp, 'baselineFrames' + YEAR + '.csv'))
    BLTs = list(csv.reader(BLTs))
    BLTs = np.asarray(BLTs)
    
    return frames, BLTs


# Load the database tables for CogCarSim logs 
#  - these must have already been exported from SQL to csv
def getCogCarSimData(inp):
    CCSsteps = np.array(np.genfromtxt(
            os.path.join(inp, 'cogcarsim2_2019_step.csv')
            , delimiter=",", skip_header = 1))
    CCSruns = open(os.path.join(inp, 'cogcarsim2_2019_run.csv'))
    CCSruns = list(csv.reader(CCSruns))
    CCSruns = np.asarray(CCSruns)
    return CCSruns, CCSsteps



# DECLARE GLOBALS
basepath = '/media/bcowley/Transcend/flowcar-backup/cogcarsim_physio_2019'
read_dir = 'ccs_recordings'
rite_dir = 'sync_physio2ccs'
sigmap = {'EDA':'E', 'EOG':['A', 'B'], 'BVP':'F'}
mapsig = {'E':'EDA', 'A':'EOG', 'B':'EOG', 'F':'BVP'}
siglist = ['EDA', 'BVP', 'EOG']
if "frames" not in globals() and "BLTs" not in globals():
    frames, BLTs = getframes(basepath)
if "CCSrun" not in globals() and "CCSstp" not in globals():
    CCSrun, CCSstp = getCogCarSimData(basepath)


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
    
    fin = open(infoplyrjsonfile)
    
    if infoplyrjsonfile.endswith('csv'):
        info = dict(csv.reader(fin))
        sync_str = "Start Time (Synced)"
        syst_str = "Start Time (System)"
        
    elif infoplyrjsonfile.endswith('json'):
        info = json.load(fin)
        sync_str = "start_time_synced_s"
        syst_str = "start_time_system_s"
    
    for k, v in info.items():
        if k == sync_str: 
            sync = v
        elif k == syst_str: 
            syst = v
    
    systime = float(syst) - float(sync)
    
    return systime


# Function to plot and save figure comparing pre- and post-interpolation data
def sanity_check(svpath, signame, t1, d1, t2, d2):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    
    ax1.plot(t1, 1.0/d1, '-')
    ax1.set_title('raw ' + signame)
    
    ax2.plot(t2, 1.0/d2, '-')
    ax2.get_yaxis().set_visible(False)
    ax2.set_title('interpolated ' + signame)
    
    if not os.path.exists(os.path.dirname(svpath)):
        os.makedirs(os.path.dirname(svpath), exist_ok = True)
        
    plt.savefig(svpath + '_interp_' + signame)
    plt.show()


# parse data path for subject, session, run indices
def get_sbj_ssn(datpth):
    
    pfx = os.sep + 'Nexus_data' + os.sep + 'session_'
    if datpth.find(pfx) == -1:
        pfx = os.sep + 'Pupil_data' + os.sep + 'session_'
        
    pfx = datpth.find(pfx) + len(pfx)
    sbj = datpth[pfx:pfx + 2]
    ssn = datpth[pfx + 3:]
    
    lg.info('Subject: ' + sbj)
    lg.info('Session: ' + ssn)
    
    return sbj, ssn


# Function parses given nexus file, 
def parse_nxs(inp, fnm, sigs):
    
    # Get world time stamp data from pupil folder
    pupilpath = os.path.join(inp.replace('Nexus', 'Pupil'), '000')
    world_times = np.load(os.path.join(pupilpath, "world_timestamps.npy"))
    world_times += getpupilsystime(os.path.join(pupilpath, "info.player.json"))

    # get the signal and timestamp data    
    chx = []
    for k, v in sigmap.items():
        if k in sigs:
            [chx.append(c) for c in v]

    hdr, dat = getsig(os.path.join(inp, fnm), chx)
    ts = dat[0,:]
    
    sbj, ssn = get_sbj_ssn(inp)

    # Get baseline segment based on BLTs
    sbj_bl = (BLTs[:,0] == sbj + '-' + ssn + '-S') | (
            BLTs[:,0] == sbj + '-' + ssn + '-E')
    bldat = []
    if sum(sbj_bl) != 2:
        lg.error('Baseline frame FU, subject {}, session {}'.format(sbj, ssn))
    else:
        frblt = np.int_(BLTs[sbj_bl, 1])
        t1 = np.argmin(abs(ts - world_times[frblt[0]]))
        t2 = np.argmin(abs(ts - world_times[frblt[1]]))
        bldat = dat[:, t1:t2]
        bldat[0, :] -= bldat[0, 0] # zero the timestamps

    # Split time and data to make interpolation of the trials readable
    dat = np.delete(dat, 0, 0)
    
    # Get correct frame/timebase for this subject, session, & run
    frssn = [len(world_times) - 1]
    tbase = [world_times[-1]]
    trls = []
    for f in np.flip(frames[:,0]):
        if (f.startswith(sbj + '-' + ssn + '-')):
            trls = [f[-2:]] + trls
            frssn = [int(frames[np.where(frames[:,0] == f),1])] + frssn
            tbase = [world_times[frssn[0]]] + tbase
    if len(frssn) == 1:
        lg.error('No trial frames for subject {}, session {}'.format(sbj, ssn))
        return hdr, trls, [], bldat, sbj, ssn
    
    outdat = [None] * (len(frssn) - 1)
    
    # Step through each run and interpolate each signal, building output data
    for i in range(len(frssn) - 1):
        lg.info('Trial {}: frame {} to frame {}'.format(
                trls[i], frssn[i], frssn[i + 1]))
        
        # Get the run times from cogcarsim database
        runidx = CCSrun[:,1] == sbj + '-0' + ssn[1] + '-' + trls[i]
        runTS = CCSstp[CCSstp[:, 1] == int(CCSrun[runidx,0]), 2].ravel()
        runTS -= runTS[0]
        # unique exception when physio data got cut off before game end
        if (inp[-6:] == "09-01A") & (trls[i] == "04"):
            runTS = runTS[0:-5081]
        lg.info("Length of run: {} seconds".format(round(runTS[-1])))
        
        # Extract biosigs and physio timestamps for trial i
        t1 = np.argmin(abs(ts - tbase[i]))
        t2 = np.argmin(abs(ts - tbase[i + 1]))
        fi_ts = ts[t1 - 128:t2] - tbase[i]
        fi_dat = dat[:, t1 - 128:t2]
        
        framei = [None] * (len(hdr) + 1)

        # extract the right physio data, and interp based on cogcarsim
        for chx in range(len(hdr)):
            lg.info("Extracting {}".format(mapsig[hdr[chx]]))
            #interpolate biosig & get final signal with cogcarsim timestamps
            interper = spi.interp1d(fi_ts, fi_dat[chx,:], kind = 1)
            try:
                framei[chx + 1] = interper(runTS)
            except:
                lg.exception("Interpolation error!")
            else:
                # plot comparison of pre- to post-interp
                sanity_check(
                        os.path.join(inp.replace(read_dir, rite_dir)
                        , 'figs', 'trial{}'.format(trls[i]))
                        , mapsig[hdr[chx]]
                        , fi_ts, fi_dat[chx,:]
                        , runTS, framei[chx + 1])

        if all([f is None for f in framei]):
            continue
        framei[0] = runTS
        outdat[i] = np.array(framei)
    
    return hdr, trls, outdat, bldat, sbj, ssn


# Function walks through folders under given path & processes any nexus files
def traversi(indir, sigs = ['EDA', 'BVP', 'EOG']):
    indir = os.path.abspath(indir)
    
    for root, subdirs, files in os.walk(indir):
        
        outdir = root.replace(read_dir, rite_dir)
        
        for f in files:
            if f.startswith('nexus') & f.endswith('.txt'):
                
                # Extract timestamps and correct columns
                try:
                    hdr, trl, dat, bl, sbj, ssn = parse_nxs(root, f, sigs)
                except:
                    lg.exception("parse_pupil error!")
                    continue

                # Manage file writing paths etc
                lg.info('Read {} from\n {}'.format(sigs, os.path.join(root, f)))
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok = True)
                lg.info('Write to\n {}'.format(outdir))
                ts = dt.datetime.fromtimestamp(
                        float(f[f.find('_') + 1:f.find('.txt')]))
                
                # Write out the Baseline
                if bl != []:
                    for ch in hdr:
                        nm = '{}_{}_baseline_Nexus{}_{}.csv'.format(sbj, ssn,
                                mapsig[ch], ts.strftime('%Y%m%d'))
                        outf = os.path.join(outdir, nm)
                        outf = open(outf, "w", newline="")
                        writer = csv.writer(outf)
                        fi_ch_bl = bl[[True] + list(hdr==ch), :]
                        writer.writerows(fi_ch_bl.T)                    

                # Write out the data
                for i in range(len(trl)):
                    if dat[i] is None:
                        continue
                    fi_dat = dat[i]
                    for ch in hdr:
                        nm = '{}_{}_trial{}_Nexus{}_{}.csv'.format(sbj, ssn, 
                              trl[i], mapsig[ch], ts.strftime('%Y%m%d'))
                        outf = os.path.join(outdir, nm)
                        outf = open(outf, "w", newline="")
                        writer = csv.writer(outf)
                        fi_ch_dat = fi_dat[[True] + list(hdr==ch), :]
                        writer.writerows(fi_ch_dat.T)

        if subdirs != []:
            for sbdr in subdirs:
                nu_ind = os.path.join(indir, sbdr)
                if os.path.exists(nu_ind):
                    lg.info('Going deeper into {}'.format(nu_ind))
                    traversi(nu_ind, sigs)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        siglist = sys.argv[1]
        if len(sys.argv) > 2:
            basepath = sys.argv[2]
            rite_dir = os.path.join(basepath, sys.argv[0])
            if len(sys.argv) > 3:
                read_dir = sys.argv[3]

# SETUP LOGGING
for handler in lg.root.handlers[:]:
    lg.root.removeHandler(handler)
lg.basicConfig(filename=os.path.join(basepath, read_dir, 'read_interp_pupil.log'), 
               filemode='w', 
               format='%(name)s - %(levelname)s - %(message)s', 
               level=lg.DEBUG)

#DEBUG
DEBUG = False
if DEBUG:
    one_ssn = '01/'
#    one_ssn = '09/Pupil_data/session_09-01/000/exports/000'
    traversi(os.path.join(basepath, read_dir, one_ssn), siglist)
else:
    traversi(os.path.join(basepath, read_dir), siglist)