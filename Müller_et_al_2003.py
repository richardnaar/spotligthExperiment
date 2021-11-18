#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.2.3),
    on November 17, 2021, at 13:43
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.2.3'
expName = 'Müller_et_al_2003'  # from the Builder filename that created this script
expInfo = {'participant': '', 'refreshRate': ['60', '120']}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\Richard Naar\\Documents\\dok\\TARU\\EEGManyLabs\\Task\\Müller_et_al_2003.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='spotlight', color=[-1,-1,-1], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='deg')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Setup eyetracking
ioDevice = ioConfig = ioSession = ioServer = eyetracker = None

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "intro"
introClock = core.Clock()
frameConst = int(expInfo['refreshRate'])/60
startFrame = 389 * frameConst
stimDur = 10.1 # (e.g. 186 + 60 frames with 60 Hz)

k = 1
boxSize = (2.5*k, 3.2*k)
symbolSize = (1.5*k, 2.5*k)

#ecc = [5.25, 10.25]
ecc = [4, 9]

xys = [(ecc[1]*-1,0), (ecc[0]*-1,0),(ecc[0],0), (ecc[1],0)]

rects = visual.ElementArrayStim(win, name = 'rects', units='deg', 
fieldPos=(0.0, 0.0), fieldSize=(24, 4), fieldShape='square', 
nElements=4, sizes= boxSize, xys=xys, 
colors=([1.0, 1.0, 1.0]) , colorSpace='rgb', opacities=1, oris=0, 
sfs=0, contrs=[1, 1,1,1], phases=0, elementTex='sqr',
elementMask=None, texRes=48, interpolate=True, 
autoLog=None, maskParams=None)

imageArray = ['stimuli/rect_ur.png','stimuli/rect_dr.png','stimuli/rect_dl.png','stimuli/rect_ul.png','stimuli/rect_target.png']
#imageArray = ['stimuli/rect_ur.jpg','stimuli/rect_dr.jpg','stimuli/rect_dl.jpg','stimuli/rect_ul.jpg','stimuli/rect_target.jpg']

randImage = [0,1,2,3,4,4,4,4]
shuffle(randImage)
nonTargetSet = imageArray[0:4]

def addData(rt, isPair, trials, nrOfEntries):
    if nrOfEntries > 0:
        thisExp.nextEntry()
    thisExp.addData('RT', rt) 
    thisExp.addData('targetPair', isPair)
    thisExp.addData('absNumOfTrials', trials)


#def loadpics(mainDir, pics, endindx, listname, units, picSize):
#    for file in range(0, endindx):
#        listname.append(visual.ImageStim(win=win, image=mainDir + '\\' +
#                                         str(pics[file]), units=units, size=picSize, name=str(pics[file])))
#stimuli = []
#loadpics(_thisDir, imageArray, len(imageArray), stimuli, 'deg', symbolSize)


# Initialize components for Routine "block_intro"
block_introClock = core.Clock()
text_block = visual.TextStim(win=win, name='text_block',
    text='',
    font='Open Sans',
    units='deg', pos=(0, 0), height=2.0, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
key_block_intro = keyboard.Keyboard()

# Initialize components for Routine "iti"
itiClock = core.Clock()
absNumOfTrials = 0
fix = visual.ShapeStim(
    win=win, name='fix', vertices='cross',units='deg', 
    size=(0.5, 0.5),
    ori=0.0, pos=(0, 0),
    lineWidth=1.0,     colorSpace='rgb',  lineColor='white', fillColor='white',
    opacity=None, depth=-1.0, interpolate=True)
ISI = clock.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='ISI')

# Initialize components for Routine "trial"
trialClock = core.Clock()
# Make keyboard object
kb = keyboard.Keyboard()
image_a = visual.ImageStim(
    win=win,
    name='image_a', units='deg', 
    image='sin', mask=None,
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
image_b = visual.ImageStim(
    win=win,
    name='image_b', units='deg', 
    image='sin', mask=None,
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)
image_c = visual.ImageStim(
    win=win,
    name='image_c', units='deg', 
    image='sin', mask=None,
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
image_d = visual.ImageStim(
    win=win,
    name='image_d', units='deg', 
    image='sin', mask=None,
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-4.0)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "intro"-------
continueRoutine = True
# update component parameters for each repeat
# keep track of which components have finished
introComponents = []
for thisComponent in introComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
introClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "intro"-------
while continueRoutine:
    # get current time
    t = introClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=introClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in introComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "intro"-------
for thisComponent in introComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# the Routine "intro" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
blocks = data.TrialHandler(nReps=1.0, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('blocks.xlsx'),
    seed=None, name='blocks')
thisExp.addLoop(blocks)  # add the loop to the experiment
thisBlock = blocks.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
if thisBlock != None:
    for paramName in thisBlock:
        exec('{} = thisBlock[paramName]'.format(paramName))

for thisBlock in blocks:
    currentLoop = blocks
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock:
            exec('{} = thisBlock[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "block_intro"-------
    continueRoutine = True
    # update component parameters for each repeat
    text_block.setText(attendCond)
    key_block_intro.keys = []
    key_block_intro.rt = []
    _key_block_intro_allKeys = []
    # keep track of which components have finished
    block_introComponents = [text_block, key_block_intro]
    for thisComponent in block_introComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    block_introClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "block_intro"-------
    while continueRoutine:
        # get current time
        t = block_introClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=block_introClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *text_block* updates
        if text_block.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_block.frameNStart = frameN  # exact frame index
            text_block.tStart = t  # local t and not account for scr refresh
            text_block.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_block, 'tStartRefresh')  # time at next scr refresh
            text_block.setAutoDraw(True)
        
        # *key_block_intro* updates
        waitOnFlip = False
        if key_block_intro.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
            # keep track of start time/frame for later
            key_block_intro.frameNStart = frameN  # exact frame index
            key_block_intro.tStart = t  # local t and not account for scr refresh
            key_block_intro.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(key_block_intro, 'tStartRefresh')  # time at next scr refresh
            key_block_intro.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(key_block_intro.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(key_block_intro.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if key_block_intro.status == STARTED and not waitOnFlip:
            theseKeys = key_block_intro.getKeys(keyList=['space'], waitRelease=False)
            _key_block_intro_allKeys.extend(theseKeys)
            if len(_key_block_intro_allKeys):
                key_block_intro.keys = _key_block_intro_allKeys[-1].name  # just the last key pressed
                key_block_intro.rt = _key_block_intro_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in block_introComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "block_intro"-------
    for thisComponent in block_introComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    blocks.addData('text_block.started', text_block.tStartRefresh)
    blocks.addData('text_block.stopped', text_block.tStopRefresh)
    # check responses
    if key_block_intro.keys in ['', [], None]:  # No response was made
        key_block_intro.keys = None
    blocks.addData('key_block_intro.keys',key_block_intro.keys)
    if key_block_intro.keys != None:  # we had a response
        blocks.addData('key_block_intro.rt', key_block_intro.rt)
    blocks.addData('key_block_intro.started', key_block_intro.tStartRefresh)
    blocks.addData('key_block_intro.stopped', key_block_intro.tStopRefresh)
    # the Routine "block_intro" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=40.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
    
    for thisTrial in trials:
        currentLoop = trials
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                exec('{} = thisTrial[paramName]'.format(paramName))
        
        # ------Prepare to start Routine "iti"-------
        continueRoutine = True
        routineTimer.add(4.000000)
        # update component parameters for each repeat
        rects.opacities = 1
        absNumOfTrials += 1
        thisExp.addData('absNumOfTrials', absNumOfTrials) 
        
        
        # keep track of which components have finished
        itiComponents = [fix, ISI]
        for thisComponent in itiComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        itiClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "iti"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = itiClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=itiClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            rects.draw()
            
            
            # *fix* updates
            if fix.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fix.frameNStart = frameN  # exact frame index
                fix.tStart = t  # local t and not account for scr refresh
                fix.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fix, 'tStartRefresh')  # time at next scr refresh
                fix.setAutoDraw(True)
            if fix.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fix.tStartRefresh + 4-frameTolerance:
                    # keep track of stop time/frame for later
                    fix.tStop = t  # not accounting for scr refresh
                    fix.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(fix, 'tStopRefresh')  # time at next scr refresh
                    fix.setAutoDraw(False)
            # *ISI* period
            if ISI.status == NOT_STARTED and t >= 0-frameTolerance:
                # keep track of start time/frame for later
                ISI.frameNStart = frameN  # exact frame index
                ISI.tStart = t  # local t and not account for scr refresh
                ISI.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(ISI, 'tStartRefresh')  # time at next scr refresh
                ISI.start(0.5)
            elif ISI.status == STARTED:  # one frame should pass before updating params and completing
                # updating other components during *ISI*
                image_a.setImage(imageArray[randImage[0]])
                image_b.setImage(imageArray[randImage[1]])
                image_c.setImage(imageArray[randImage[2]])
                image_d.setImage(imageArray[randImage[3]])
                # component updates done
                ISI.complete()  # finish the static period
                ISI.tStop = ISI.tStart + 0.5  # record stop time
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in itiComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "iti"-------
        for thisComponent in itiComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        trials.addData('fix.started', fix.tStartRefresh)
        trials.addData('fix.stopped', fix.tStopRefresh)
        trials.addData('ISI.started', ISI.tStart)
        trials.addData('ISI.stopped', ISI.tStop)
        
        # ------Prepare to start Routine "trial"-------
        continueRoutine = True
        # update component parameters for each repeat
        counter = 0
        switchTime = 0
        frameCount = 0
        responseGiven = False
        targetPair = False
        nrOfEntries = 0
        cond = [pos2attLeft-1, pos2attRight-1]
        
        imList = [image_a, image_b, image_c, image_d]
        
        #image_a.opacity=1
        #image_b.opacity=1
        #image_c.opacity=1
        #image_d.opacity=1
        
        
        image_a.setContrast(1.0)
        image_a.setPos([xys[0]])
        image_a.setSize(symbolSize)
        image_b.setPos([xys[1]])
        image_b.setSize(symbolSize)
        image_c.setPos([xys[2]])
        image_c.setSize(symbolSize)
        image_d.setPos([xys[3]])
        image_d.setSize(symbolSize)
        # keep track of which components have finished
        trialComponents = [image_a, image_b, image_c, image_d]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        trialClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
        frameN = -1
        
        # -------Run Routine "trial"-------
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            tThisFlip = win.getFutureFlipTime(clock=trialClock)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            rects.draw()
            fix.draw()
            
            if t < 9.1:
                opacities = [0,0,0,0]
                frameNow = startFrame + frameN
            
                if frameNow % (4*frameConst) == 0:
                    opacities[0] = 1
            
                if frameNow % (7*frameConst) == 0:
                    opacities[1] = 1
            
                if frameNow % (3*frameConst) == 0:
                    opacities[2] = 1
            
                if frameNow % (5*frameConst) == 0:
                    opacities[3] = 1
            
            #    rects.opacities = opacities
            
                if frameNow % (6*11*frameConst) == 0:
                    shuffle(randImage)
                    imCount = 0
                    for im in imList:
                        im.image = imageArray[randImage[imCount]]
                        imCount += 1
                    if 'target' in imageArray[randImage[cond[0]]] and 'target' in imageArray[randImage[cond[1]]]:
                        if frameN-frameCount > 2*54 or not nrOfEntries:
                            if targetPair and not responseGiven:
                                addData('noResp', targetPair, absNumOfTrials, nrOfEntries)
                                nrOfEntries += 1
                            switchTime = t
                            targetPair = True
                            frameCount = frameN
                            responseGiven = False
                        else:
                            imList[cond[randint(0,1)]].image = nonTargetSet[randint(0,4)]
            
                image_a.draw()
                image_b.draw()
                image_c.draw()
                image_d.draw()
            
                keys = kb.getKeys()
            
                if keys:
                    if not targetPair and 'space' in keys[-1].name:
                        addData(t-switchTime, targetPair, absNumOfTrials, nrOfEntries)
                        nrOfEntries += 1
                        responseGiven = True
                    elif 'escape' in keys[-1].name:
                        core.quit()
            
            
                if keys and not responseGiven:
                    responseGiven = True
                    if t-switchTime > 0.25 and t-switchTime <= 1:
                        if 'space' in keys[-1].name:
                            addData(t-switchTime, targetPair, absNumOfTrials, nrOfEntries)
                            nrOfEntries += 1
                            responseGiven = True
                    elif ('space' in keys[-1].name and targetPair) or ('space' in keys[-1].name and t-switchTime <= 0.25):
                            addData(t-switchTime, targetPair, absNumOfTrials, nrOfEntries)
                            responseGiven = True
            #                targetPair = False
                            nrOfEntries += 1
            #else:
            #    image_a.opacity=0
            #    image_b.opacity=0
            #    image_c.opacity=0
            #    image_d.opacity=0
            #    rects.opacities = 1
            
            
            
            # *image_a* updates
            if image_a.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image_a.frameNStart = frameN  # exact frame index
                image_a.tStart = t  # local t and not account for scr refresh
                image_a.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image_a, 'tStartRefresh')  # time at next scr refresh
                image_a.setAutoDraw(True)
            if image_a.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image_a.tStartRefresh + stimDur-frameTolerance:
                    # keep track of stop time/frame for later
                    image_a.tStop = t  # not accounting for scr refresh
                    image_a.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(image_a, 'tStopRefresh')  # time at next scr refresh
                    image_a.setAutoDraw(False)
            
            # *image_b* updates
            if image_b.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image_b.frameNStart = frameN  # exact frame index
                image_b.tStart = t  # local t and not account for scr refresh
                image_b.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image_b, 'tStartRefresh')  # time at next scr refresh
                image_b.setAutoDraw(True)
            if image_b.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image_b.tStartRefresh + stimDur-frameTolerance:
                    # keep track of stop time/frame for later
                    image_b.tStop = t  # not accounting for scr refresh
                    image_b.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(image_b, 'tStopRefresh')  # time at next scr refresh
                    image_b.setAutoDraw(False)
            
            # *image_c* updates
            if image_c.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image_c.frameNStart = frameN  # exact frame index
                image_c.tStart = t  # local t and not account for scr refresh
                image_c.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image_c, 'tStartRefresh')  # time at next scr refresh
                image_c.setAutoDraw(True)
            if image_c.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image_c.tStartRefresh + stimDur-frameTolerance:
                    # keep track of stop time/frame for later
                    image_c.tStop = t  # not accounting for scr refresh
                    image_c.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(image_c, 'tStopRefresh')  # time at next scr refresh
                    image_c.setAutoDraw(False)
            
            # *image_d* updates
            if image_d.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image_d.frameNStart = frameN  # exact frame index
                image_d.tStart = t  # local t and not account for scr refresh
                image_d.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image_d, 'tStartRefresh')  # time at next scr refresh
                image_d.setAutoDraw(True)
            if image_d.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image_d.tStartRefresh + stimDur-frameTolerance:
                    # keep track of stop time/frame for later
                    image_d.tStop = t  # not accounting for scr refresh
                    image_d.frameNStop = frameN  # exact frame index
                    win.timeOnFlip(image_d, 'tStopRefresh')  # time at next scr refresh
                    image_d.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # the Routine "trial" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 40.0 repeats of 'trials'
    
    thisExp.nextEntry()
    
# completed 1.0 repeats of 'blocks'


# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
