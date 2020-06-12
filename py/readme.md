
# Timestamp syncing

## 2017 procedure

First, trial start frames of Pupil videos were estimated and logged in the Frames.csv file (instructions below).

__Searching for frames in pupil__:

1. Open the Pupil Player on the Linux lab laptop
    - Go to: ~/pupil/pupil\_src/player/
    - Do it: $ python3 main.py
2. grey Pupil player window opens → drag data folder into window (e.g. 000)
3. pause video around the time the game is started, move between frames with left/right arrow keys
4. target frame is the first one where the game text is not visible
5. frame is the number next to the play/pause button
6. write down run (e.g. 01-1-1), frame number, and pupil file


The resulting Frames.csv file was used in the interpolate\_physiology.py script to set trial start times (unix timestamps) for physiology data. The physiology data was then interpolated to match cogcarsim timestamps.

 

The script requires a bit of "manual work" (i.e. does not loop over trials), due to our limited programming skills at the time, and the  need for checking each recording individually anyway. It is assumed that there is one recording file per trial (and that files are named after the trials: 01-1-1.jsons, and so on) but that was always not the case, so it might cause an error. However it should be seen from filenames (e.g. participant 5 session 5 had all trials in one recording if I remember correctly).

 

Also, some recordings had missing end/start segments so they had to be handled differently. physiology\_notes.ods has the missing trial/segment info.  
