import psychopy
from psychopy import locale_setup, gui, visual, core, data, event, logging, monitors, sound
import os

from numpy.random import random, randint, shuffle, seed

import numpy as np

expName = os.path.basename(__file__)
expInfo = {'Participant': '001'}

# dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
# if dlg.OK == False:
#     core.quit()  # user pressed cancel

expInfo['date'] = data.getDateStr()  # add a simple timestamp

# durations
trialDuration = 3.060  # 3060 ms


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
gapPos = [(rectanglePos[0]+horiz/2, 0.56), (rectanglePos[1]+horiz/2, -0.56),
          (rectanglePos[2]-horiz/2, 0.56), (rectanglePos[3], 0), (rectanglePos[3]-horiz, 0)]
posDic = {'1': [gapPos[0], 0], '2': [gapPos[0], 0], '3': [
    gapPos[0], 0], '4': [gapPos[0], 0], '5': [gapPos[0], 90]}

# rectangle with a right upper gap
box_1 = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(rectanglePos[0], 0),
    lineWidth=16, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

gap_1 = visual.Rect(
    win=win, units='deg',
    width=(1, 1)[0], height=(vert/2, vert/2)[1],
    ori=0, pos=(rectanglePos[0]+horiz/2, 0.56),
    lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# rectangle with a gap in right and down
box_2 = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(rectanglePos[1], 0),
    lineWidth=16, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

gap_2 = visual.Rect(
    win=win, units='deg',
    width=(1, 1)[0], height=(vert/2, vert/2)[1],
    ori=0, pos=(rectanglePos[1]+horiz/2, -0.56),
    lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# rectangle with a left upper gap
box_3 = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(rectanglePos[2], 0),
    lineWidth=16, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

gap_3 = visual.Rect(
    win=win, units='deg',
    width=(1, 1)[0], height=(vert/2, vert/2)[1],
    ori=0, pos=(rectanglePos[2]-horiz/2, 0.56),
    lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# rectangle with a gap in left and down
box_4 = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(rectanglePos[3], 0),
    lineWidth=16, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

gap_4 = visual.Rect(
    win=win, units='deg',
    width=(1, 1)[0], height=(vert/2, vert/2)[1],
    ori=0, pos=(rectanglePos[3]-horiz/2, -0.56),
    lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

box_5 = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(rectanglePos[3], 0),
    lineWidth=16, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

# gap_5 = visual.Rect(
#     win=win, units='deg',
#     width=(1, 1)[0], height=(vert/2, vert/2)[1],
#     ori=90, pos=(rectanglePos[3]-horiz, 0),
#     lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
#     fillColor=[1, 1, 1], fillColorSpace='rgb',
#     opacity=1, depth=0.0, interpolate=True)

fixation = visual.ShapeStim(
    win=win, name='fixation', vertices='cross',
    size=(.8, .8),
    ori=0, pos=(0, 0),
    lineWidth=0, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[1, 1, 1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

runExperiment = True

trialNumber = 0
while runExperiment:
    trialNumber += 1

    box_1.draw(), gap_1.draw()
    box_2.draw(), gap_2.draw()

    box_3.draw(), gap_3.draw()
    box_4.draw(), gap_4.draw()

    fixation.draw()

    win.flip()

core.wait(3)

win.close(), core.quit()

# The symbols were presented simultaneously at the 4 locations for 11 frames (180.95 ms)
# followed immediately by the next symbol array.
# Thus, the symbol onsets and offsets did
# not occur in synchrony with the background flickering rectangles that drove the SSVEP.
# Target symbols occurred equally often at all stimulus positions, following randomized
# sequences.

# read xls table in
xls_file = pd.ExcelFile(fileName + '.xlsx')
table = xls_file.parse()  # pars to table

# sort table by emo and picset
table = table.sort_values(['emo', 'picset']).reset_index(drop=True)

# shuffle a whole table
shuffledTable = table.sample(frac=1).reset_index(drop=True)

# randomize two columns only
table[['cond', 'presentVAS']] = table[['cond', 'presentVAS']].sample(
    frac=1).reset_index(drop=True)

np.unique(table['emo'])

# shuffle a list
somelist = list(range(1, 10))
np.random.shuffle(somelist)

# repeat the list 4 times
condlist = somelist*4

# 2 by 10 arrey
np.zeros((2, 10))

# 2 random integers without replacement from 0-9
np.random.choice(10, 2, replace=False)
