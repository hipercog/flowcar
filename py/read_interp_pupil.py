#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:27:26 2022

@author: bcowley
"""

import os, json, csv, sys
import numpy as np
import pandas as pd
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
siglist = ['confidence', 'norm_pos_x', 'norm_pos_y', 'diameter']
sigix = [0, 1, 2, 3]
if "frames" not in globals() and "BLTs" not in globals():
    frames, BLTs = getframes(basepath)
if "CCSrun" not in globals() and "CCSstp" not in globals():
    CCSrun, CCSstp = getCogCarSimData(basepath)



## Get signal out of JSON file
#def getsig(ppl_pos_file, chans):
#    ppldata = []
#    
#    #    create default data in case there are reading problems
#    ppl_ts_i = {'ts':math.nan}
#    ppl_sg_i = dict((k,math.nan) for k in chans)
#    
#    with open(ppl_pos_file, newline='') as f:
#        i = 0
#        ppl_rdr = csv.reader(f, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
#        for line in ppl_rdr:
#            i += 1
#            try:
#                ppl_ts_i = line[0]
#                ppl_sg_i = line[1]
#            except:
#                lg.exception("csv read error line{}: {}".format(str(i), line))
#            else:
#                fields = [ppl_ts_i[1]]
#                for c in chans:
#                    fields.append(ppl_sg_i[c])
#                ppldata.append(fields)
#            
#        f.close()
#        
#    pplarr = np.array(ppldata).T
#    idx = np.divide([len(np.unique(r)) * 100 for r in pplarr], len(pplarr.T)) > 5
#        
#    return np.array(chans)[idx[1:]], pplarr[idx,:]


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
    
    pfx = os.sep + 'Pupil_data' + os.sep + 'session_'
    if datpth.find(pfx) == -1:
        pfx = os.sep + 'Nexus_data' + os.sep + 'session_'
        
    pfx = datpth.find(pfx) + len(pfx)
    sbj = datpth[pfx:pfx + 2]
    ssn = datpth[pfx + 3:pfx + 5]
    
    lg.info('Subject: ' + sbj)
    lg.info('Session: ' + ssn)
    
    return sbj, ssn


# Function parses given pupil_positions file, 
def parse_ppl(inp, fnm, sigs):
    
    # Get world time stamp data from pupil folder
    pupilpath = inp.replace(os.path.join('exports', '000'), '')
    world_times = np.load(os.path.join(pupilpath, "world_timestamps.npy"))
    pupilpath = os.path.join(inp, fnm)

    # get the signal and timestamp data
    ppldata = pd.read_csv(pupilpath)
    hdr = np.array([siglist[i] for i in sigs])
    dat = ppldata.loc[:,np.isin(ppldata.columns, hdr)].to_numpy().T
    ts = ppldata.iloc[:,0].to_numpy().T
    
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
            lg.info("Extracting {}".format(hdr[chx]))
            #interpolate biosig & get final signal with cogcarsim timestamps
            interper = spi.interp1d(fi_ts, fi_dat[chx,:], kind = 1)
            try:
                framei[chx + 1] = interper(runTS)
            except:
                lg.exception("Interpolation error!")
            else:
                # plot comparison of pre- to post-interp
                sanity_check(
                        os.path.join(inp.replace(read_dir, rite_dir).replace(
                                os.path.join('000', 'exports', '000'), '')
                                , 'figs', 'trial{}'.format(trls[i]))
                        , hdr[chx]
                        , fi_ts, fi_dat[chx,:]
                        , runTS, framei[chx + 1])

        if all([f is None for f in framei]):
            continue
        framei[0] = runTS
        outdat[i] = np.array(framei)
    
    return hdr, trls, outdat, bldat, sbj, ssn

def get_pupil_timestamp(inpath):
    inp = inpath.replace('Pupil', 'Nexus').replace(
                        os.path.join('000', 'exports', '000'), '')
    for rt, sb, fyles in os.walk(inp):
        for fy in fyles:
            if fy.startswith('pupil_info_'):
                ts = dt.datetime.fromtimestamp(
                    float(fy[fy.find('o_') + 2:fy.find('.txt')]))
                return ts


# Function walks through folders under given path & processes any pupil files
def traversi(inpath, sigs):
    inpath = os.path.abspath(inpath)
    
    for root, subdirs, files in os.walk(inpath):
                
        for f in files:
            if f == 'pupil_positions.csv':
                outdir = root.replace(os.path.join('000', 'exports', '000'), ''
                                    ).replace(read_dir, rite_dir)
                ts = get_pupil_timestamp(root)
                if ts is None:
                    print('No timestamp in {}'.format(root))
                
                # Extract timestamps and correct columns
#                FIXME ADD TIMESTAMPS TO OUTPUT SO INDIVIDUAL CHANNEL FILES CAN HAVE
                try:
                    hdr, trl, dat, bl, sbj, ssn = parse_ppl(root, f, sigs)
                    if dat == []:
                        continue
                except:
                    lg.exception("parse_pupil error!")
                    continue

                # Manage file writing paths etc
                lg.info('Read {} from\n{}'.format(hdr, os.path.join(root, f)))
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok = True)
                lg.info('Write to\n{}'.format(outdir))
                
                # Write out the Baseline
                if bl != []:
                    for ch in hdr:
                        nm = '{}_{}_baseline_Pupil{}_{}.csv'.format(sbj, ssn,
                                ch, ts.strftime('%Y%m%d'))
                        outf = os.path.join(outdir, nm)
                        outf = open(outf, "w", newline="")
                        writer = csv.writer(outf)
                        fi_ch_bl = bl[np.isin(hdr, ch), :]
                        writer.writerows(fi_ch_bl.T)

                # Write out the data
                for i in range(len(trl)):
                    if dat[i] is None:
                        continue
                    fi_dat = dat[i]
                    for ch in hdr:
                        nm = '{}_{}_trial{}_Pupil{}_{}.csv'.format(sbj, ssn, 
                              trl[i], ch, ts.strftime('%Y%m%d'))
                        outf = os.path.join(outdir, nm)
                        outf = open(outf, "w", newline="")
                        writer = csv.writer(outf)
                        fi_ch_dat = fi_dat[[True] + list(np.isin(hdr, ch)), :]
                        writer.writerows(fi_ch_dat.T)
#                f = [] # Clear f because it might be the last file

        if subdirs != []:
            for sbdr in subdirs:
                nu_ind = os.path.join(root, sbdr)
                if os.path.exists(nu_ind):
                    lg.info('Going deeper into {}'.format(nu_ind))
                    traversi(nu_ind, sigs)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sigix = sys.argv[1]
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
DEBUG = True
if DEBUG:
    one_ssn = '09/'
#    one_ssn = '09/Pupil_data/session_09-01/000/exports/000'
    traversi(os.path.join(basepath, read_dir, one_ssn), sigix)
else:
    traversi(os.path.join(basepath, read_dir), sigix)