import msgpack
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
import scipy.interpolate

rec_dir = "/home/tru/pupil/recordings/2017_10_20/000/"

world_timestamps = np.load(os.path.join(rec_dir, "world_timestamps_unix.npy"))
frame = 200
timebase = world_timestamps[frame]
print(timebase)

f = open(os.path.join(rec_dir, "pupil_data"), 'rb')
data = msgpack.unpack(f)
gaze_positions = data['gaze_positions']
gaze_norm = [[g['base_data'][0]['unix_ts'] - timebase, g['norm_pos'][0], g['norm_pos'][1]] for g in gaze_positions]
gaze_norm = np.array(gaze_norm)
np.savetxt("data.csv", gaze_norm, delimiter = ",")

print(len(gaze_norm[:,0]))

data = genfromtxt("/home/tru/pupil/recordings/2017_10_20/000/step.csv", delimiter=",", skip_header = 1)

#gaze_norm = gaze_norm[(gaze_norm[:,0] > cogcarsimdata[0,2]) & ]


cogcarsimdata = np.array(data)

#print cogcarsimdata[:,2]

f = scipy.interpolate.interp1d(gaze_norm[:,0], gaze_norm[:,1], kind = 1)

final_signal = f(cogcarsimdata[:,2])

print(len(cogcarsimdata[:,0]))

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(gaze_norm[:,0], gaze_norm[:,1], 'o')
ax1.set_title('raw horz gaze')

#plt.figure()

ax2.plot(cogcarsimdata[:,2], final_signal, 'o')
ax2.set_title('cogcarsim')

plt.show()


