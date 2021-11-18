# 22/09/2020
# each location flickering at a different rate
# 15.2 Hz (position 1), 8.7 Hz (position 2), 20.3 Hz (position 3) and 12.2 Hz (position 4).
# Symbol presentations occurred in synchrony at the four locations with fixed durations of 181 ms.
# import psychopy
# locale_setup, gui, logging, monitors, sound
from psychopy import visual, core, data, event
import os

# from numpy.random import random, randint, shuffle, seed

import numpy as np

expName = os.path.basename(__file__)
expInfo = {'Participant': '001'}

# dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
# if dlg.OK == False:
#     core.quit()  # user pressed cancel

# get the current directory
dirpath = os.getcwd()
stimDir = dirpath + '\\stimuli'

expInfo['date'] = data.getDateStr()  # add a simple timestamp

# durations
trialDuration = 3.060  # 3060 ms

clock = core.Clock()
# 12 blocks of 40 trials each, resulting in 120 trials per experimental condition.
# The responding hand was changed halfway through the experiment, and the sequence of hand usage was counterbalanced across subjects.

# cycle durarions of 2,7,3,5 frames (on one frame and followed apropriate number of off frames)

monSettings = {'size': (800, 600), 'fullscr': False}

win = visual.Window(
    size=monSettings['size'], fullscr=monSettings['fullscr'], screen=0, color='black',
    blendMode='avg', useFBO=False, monitor='testMonitor',
    units='deg', waitBlanking=True)

# Five different symbols were used, with one symbol serving as the target.
frex = [15.20, 8.69, 20.27, 12.16]
horiz = 2.5
vert = 3.2

rectanglePos = [-9, -4, 4, 9]

rects = list(filter(lambda x: x.endswith('.jpg'), os.listdir(stimDir)))

#listname = images[indx]
# current_image.draw()

fixation = visual.ShapeStim(
    win=win, name='fixation', vertices='cross',
    size=(1, 1),
    ori=0, pos=(0, 0),
    lineWidth=1, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[1, 1, 1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)


def draw_fix(win, fixation, duration):
    fixStartTime = clock.getTime()
    time = clock.getTime() - fixStartTime
    fixPresented = False
    while (time) < duration:
        if not event.getKeys('q'):

            fixation.draw()
            win.flip()
            if not fixPresented:
                fixPresented = True
        else:
            core.quit()
        time = clock.getTime() - fixStartTime


def loadpics(picture_directory, pics, endindx, listname, units, picSize):
    for file in range(0, endindx):
        listname.append(visual.ImageStim(win=win, image=picture_directory + '\\' +
                                         str(pics[file]), units=units, size=picSize, name=str(pics[file])))


def randomStim(randStim):
    randStim = list(np.random.choice(5, 4, replace=True))


stimuli = []
loadpics(stimDir, rects, len(rects), stimuli, 'deg', (horiz, vert))
#listname = images[indx]
# current_image.draw()

# On separate blocks of trials, subjects were instructed
# verbally to attend to either the two left field positions (1 + 2), the
# two right field positions (3 + 4), or to two separated positions
# (1 + 3) or (2 + 4).
# Simultaneous target symbols occurred unpre dictably at the two attended locations (0–3 times per trial), as well as
# at the other locations.


runExperiment = True

trialNumber, frameN = 0, 0
randStim = list(np.random.choice(5, 4, replace=True))
expTime = clock.getTime()
while runExperiment:
    trialNumber += 1
    frameN += 1

    if not event.getKeys('q'):
        runExperiment = True
    else:
        core.quit()

    # if frameN % 11 == 0:
    #     randStim = list(np.random.choice(5, 4, replace=True))
    # randStim = list(np.random.choice(5, 4, replace=True))

    if frameN % 11 == 0:

        randomStim(randStim)
        fixation.draw()

        stimuli[randStim[0]].pos = (rectanglePos[0], 0)
        stimuli[randStim[0]].draw()

        stimuli[randStim[1]].pos = (rectanglePos[1], 0)
        stimuli[randStim[1]].draw()

        stimuli[randStim[2]].pos = (rectanglePos[2], 0)
        stimuli[randStim[2]].draw()

        stimuli[randStim[3]].pos = (rectanglePos[3], 0)
        stimuli[randStim[3]].draw()

        win.flip()

    if clock.getTime() - expTime > trialDuration:
        fixation.draw()
        win.flip()
        core.wait(1.5)
        expTime = clock.getTime()

    # draw_fix(win, fixation, 1)


win.close(), core.quit()

# The symbols were presented simultaneously at the 4 locations for 11 frames (180.95 ms)
# followed immediately by the next symbol array.
# Thus, the symbol onsets and offsets did
# not occur in synchrony with the background flickering rectangles that drove the SSVEP.
# Target symbols occurred equally often at all stimulus positions, following randomized
# sequences.
#from psychopy import visual
