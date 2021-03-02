import numpy as np
from numpy import genfromtxt
import os
import sys

#ohjeita: tarkista, etta tama tiedosto ja ts_session_6.csv ovat samassa kansiossa (Downloads)
#Downloads-kansiossa klikkaa oikealla, run in terminal -> 'python ts_frame.py'
#framet tallentuu uuteen frames_session_6.csv-tiedostoon 4. sarakkeeseen
#timestampit on 3-minuuttisten alkuaikoja. jos tarvitaan loppumisaikojen framet, 
# tulee ts_session_6.csv-tiedoston timestamp-kohtiin lisata +240
#muista tarkistaa myös Transcendin kansiopolku (rec_dir)

rec_dir = "/media/tru/Transcend/FLOW/participant_data/"



#go through each row of timestamp file 
#get frames within a window (-0.025...0.025s) around the timestamp, select first of them
framelist = []
timestamps = genfromtxt("/home/tru/Downloads/ts_session_6.csv", delimiter=",", dtype=('S15'),skip_header=1)
timestamps = np.array(timestamps)
for i in range(0,9):
	participant = timestamps[i][0]
	timestamp = float(timestamps[i][1])
	folder = timestamps[i][2]
	world_timestamps = np.load(os.path.join(rec_dir, participant, "/session_6/eye_tracking/", folder, "/world_timestamps_unix.npy"))
	frames = np.where((world_timestamps<(timestamp+0.025)) & (world_timestamps >(timestamp-0.025)))
	frame = np.concatenate(frames)[0]
	framelist.append([frame])

framelist = np.array(framelist)

#save frames to new file (last column)
ts_frames = np.hstack((timestamps.reshape(9,-1), framelist.reshape(9,-1)))
np.savetxt("frames_session_6.csv", ts_frames, delimiter = ',',fmt="%s")
