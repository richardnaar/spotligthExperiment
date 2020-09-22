import psychopy
from psychopy import locale_setup, gui, visual, core, data, event, logging, monitors, sound
import os

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
rectanglePos = [-9, -4, 4, 9]

monSettings = {'size': (800, 600), 'fullscr': False}

win = visual.Window(
    size=monSettings['size'], fullscr=monSettings['fullscr'], screen=0, color='black',
    blendMode='avg', useFBO=False, monitor='testMonitor',
    units='deg', waitBlanking=True)

# Five different symbols were used, with one symbol serving as the target.
frex = [15.20, 8.69, 20.27, 12.16]
horiz = 2.5
vert = 3.2
box = visual.Rect(
    win=win, units='deg',
    width=(horiz, horiz)[0], height=(vert, vert)[1],
    ori=0, pos=(0, 0),
    lineWidth=10, lineColor=[1, 1, 1], lineColorSpace='rgb',
    fillColor=[-1, -1, -1], fillColorSpace='rgb',
    opacity=1, depth=0.0, interpolate=True)

box.draw()
win.flip()
core.wait(3)

win.close(), core.quit()

# The symbols were presented simultaneously at the 4 locations for 11 frames (180.95 ms)
# followed immediately by the next symbol array.
# Thus, the symbol onsets and offsets did
# not occur in synchrony with the background flickering rectangles that drove the SSVEP.
# Target symbols occurred equally often at all stimulus positions, following randomized
# sequences.


# BASIC EDITING
# Alt+ ↑ / ↓ # Move line up/down
# Shift+Alt + ↓ / ↑ # Copy line up/down
# Ctrl+Shift+K # Delete line
# Ctrl+Shift+\ # Jump to matching br
# Crtl+/ # comment/un-comment

# FOLDING
# Ctrl+Shift+[ # Fold (collapse) region
# Ctrl+Shift+] # Unfold (uncollapse) region
# Ctrl+K Ctrl+[  # Fold (collapse) all subregions
# Ctrl+K Ctrl+] # Unfold (uncollapse) all subregions
