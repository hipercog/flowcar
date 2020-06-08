import glob
import csv
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import sys
from numpy import genfromtxt
import scipy.interpolate
import csv

frames = open("/home/tru/pupil_interpolation/Frames.csv")
frames = list(csv.reader(frames))
#change input here: row of csv file
row = 81
  
id = frames[row-1][0][0:2] 
session_run = frames[row-1][0][3:6]
frame = int(frames[row-1][1]) 
current_run = int(frames[row-1][2])  

print('index: ' + id + "-" + session_run)
print 'frame: ' + str(frame)
print('run: ' + str(current_run))

rec_dir = "/home/tru/pupil_interpolation/" + id

#get timebase (start time of trial) from Pupil frame
world_timestamps = np.load(os.path.join(rec_dir, session_run + "/world_timestamps_unix.npy"))
timebase = world_timestamps[frame]
print('timebase: ' + str(timebase))

#get physiological data
eda = []
for line in open(os.path.join(rec_dir, id + "-" + session_run +".jsons")):
    hdr, o = json.loads(line)
    eda.append([hdr['ts'], o['E']])

ts, eda = np.array(eda).T

ts -= timebase


bvp = []
for line in open(os.path.join(rec_dir, id + "-" + session_run +".jsons")):
    hdr, o = json.loads(line)
    bvp.append([hdr['ts'], o['F']])

ts, bvp = np.array(bvp).T
print ts
ts -= timebase


#get cogcarsim data. step is from cogcarsim2.db
data = genfromtxt("/home/tru/pupil_interpolation/step.csv", delimiter=",", skip_header = 1)
cogcarsimdata = np.array(data)
#for this run:
rundata = cogcarsimdata[np.where(cogcarsimdata[:,0]==current_run),1].ravel()
rundata -= rundata[0]
print "length of run: " + str(rundata[len(rundata)-1])


#interpolate EDA
f_eda = scipy.interpolate.interp1d(ts, eda, kind = 1)
#final EDA signal with cogcarsim timestamps
final_eda = f_eda(rundata)
#plot raw & interpolated EDA
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(ts, 1.0/eda, '-')
#ax1.set_xlim(0,ts[len(ts)-1])
ax1.set_title('raw eda')

ax2.plot(rundata, 1.0/final_eda, '-')
#ax2.set_xlim(0, ts[len(ts)-1])
ax2.set_title('interpolated eda')

plt.show()

#interpolate BVP
f_bvp = scipy.interpolate.interp1d(ts, bvp, kind = 1)
#final bvp signal with cogcarsim timestamps
final_bvp = f_bvp(rundata)
#plot raw & interpolated bvp
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(ts, bvp, '-')
#ax1.set_xlim(0,ts[len(ts)-1])
ax1.set_title('raw bvp')

ax2.plot(rundata, final_bvp, '-')
#ax2.set_xlim(0, ts[len(ts)-1])
ax2.set_title('interpolated bvp')

plt.show()
