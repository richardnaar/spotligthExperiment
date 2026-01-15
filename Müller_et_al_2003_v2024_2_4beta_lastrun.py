#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2025.2.4),
    on jaanuar 15, 2026, at 17:26
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

import psychopy
psychopy.useVersion('2025.2.4')


# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware
from psychopy.tools import environmenttools
from psychopy.constants import (
    NOT_STARTED, STARTED, PLAYING, PAUSED, STOPPED, STOPPING, FINISHED, PRESSED, 
    RELEASED, FOREVER, priority
)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard

# --- Setup global variables (available in all functions) ---
# create a device manager to handle hardware (keyboards, mice, mirophones, speakers, etc.)
deviceManager = hardware.DeviceManager()
# ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
# store info about the experiment session
psychopyVersion = '2025.2.4'
expName = 'Müller_et_al_2003_v2024_2_4beta'  # from the Builder filename that created this script
expVersion = ''
# a list of functions to run when the experiment ends (starts off blank)
runAtExit = []
# information about this experiment
expInfo = {
    'participant': f"{randint(0, 999999):06.0f}",
    'refreshRate': ['60','120'],
    'translation': ['ENG','EST'],
    'EEG': ['0','1'],
    'port address': '0x3FF8',
    'testRun': ['0','1'],
    'volume': '0.25',
    'date|hid': data.getDateStr(),
    'expName|hid': expName,
    'expVersion|hid': expVersion,
    'psychopyVersion|hid': psychopyVersion,
}

# --- Define some variables which will change depending on pilot mode ---
'''
To run in pilot mode, either use the run/pilot toggle in Builder, Coder and Runner, 
or run the experiment with `--pilot` as an argument. To change what pilot 
#mode does, check out the 'Pilot mode' tab in preferences.
'''
# work out from system args whether we are running in pilot mode
PILOTING = core.setPilotModeFromArgs()
# start off with values from experiment settings
_fullScr = True
_winSize = [1920, 1080]
# if in pilot mode, apply overrides according to preferences
if PILOTING:
    # force windowed mode
    if prefs.piloting['forceWindowed']:
        _fullScr = False
        # set window size
        _winSize = prefs.piloting['forcedWindowSize']
    # replace default participant ID
    if prefs.piloting['replaceParticipantID']:
        expInfo['participant'] = 'pilot'

def showExpInfoDlg(expInfo):
    """
    Show participant info dialog.
    Parameters
    ==========
    expInfo : dict
        Information about this experiment.
    
    Returns
    ==========
    dict
        Information about this experiment.
    """
    # show participant info dialog
    dlg = gui.DlgFromDict(
        dictionary=expInfo, sortKeys=False, title=expName, alwaysOnTop=True
    )
    if dlg.OK == False:
        core.quit()  # user pressed cancel
    # return expInfo
    return expInfo


def setupData(expInfo, dataDir=None):
    """
    Make an ExperimentHandler to handle trials and saving.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    dataDir : Path, str or None
        Folder to save the data to, leave as None to create a folder in the current directory.    
    Returns
    ==========
    psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    # remove dialog-specific syntax from expInfo
    for key, val in expInfo.copy().items():
        newKey, _ = data.utils.parsePipeSyntax(key)
        expInfo[newKey] = expInfo.pop(key)
    
    # data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
    if dataDir is None:
        dataDir = _thisDir
    filename = u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])
    # make sure filename is relative to dataDir
    if os.path.isabs(filename):
        dataDir = os.path.commonprefix([dataDir, filename])
        filename = os.path.relpath(filename, dataDir)
    
    # an ExperimentHandler isn't essential but helps with data saving
    thisExp = data.ExperimentHandler(
        name=expName, version=expVersion,
        extraInfo=expInfo, runtimeInfo=None,
        originPath='C:\\Users\\pcadmin\\Documents\\dok\\TARU\\PROJEKTID\\EEGManyLabs\\Task\\spotlight_replication\\Müller_et_al_2003_v2024_2_4beta_lastrun.py',
        savePickle=True, saveWideText=True,
        dataFileName=dataDir + os.sep + filename, sortColumns='time'
    )
    thisExp.setPriority('thisRow.t', priority.CRITICAL)
    thisExp.setPriority('expName', priority.LOW)
    # return experiment handler
    return thisExp


def setupLogging(filename):
    """
    Setup a log file and tell it what level to log at.
    
    Parameters
    ==========
    filename : str or pathlib.Path
        Filename to save log file and data files as, doesn't need an extension.
    
    Returns
    ==========
    psychopy.logging.LogFile
        Text stream to receive inputs from the logging system.
    """
    # set how much information should be printed to the console / app
    if PILOTING:
        logging.console.setLevel(
            prefs.piloting['pilotConsoleLoggingLevel']
        )
    else:
        logging.console.setLevel('warning')
    # save a log file for detail verbose info
    logFile = logging.LogFile(filename+'.log')
    if PILOTING:
        logFile.setLevel(
            prefs.piloting['pilotLoggingLevel']
        )
    else:
        logFile.setLevel(
            logging.getLevel('info')
        )
    
    return logFile


def setupWindow(expInfo=None, win=None):
    """
    Setup the Window
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    win : psychopy.visual.Window
        Window to setup - leave as None to create a new window.
    
    Returns
    ==========
    psychopy.visual.Window
        Window in which to run this experiment.
    """
    if PILOTING:
        logging.debug('Fullscreen settings ignored as running in pilot mode.')
    
    if win is None:
        # if not given a window to setup, make one
        win = visual.Window(
            size=_winSize, fullscr=_fullScr, screen=0,
            winType='pyglet', allowGUI=False, allowStencil=False,
            monitor='spotlight', color=[-1,-1,-1], colorSpace='rgb',
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=True,
            units='deg',
            checkTiming=False  # we're going to do this ourselves in a moment
        )
    else:
        # if we have a window, just set the attributes which are safe to set
        win.color = [-1,-1,-1]
        win.colorSpace = 'rgb'
        win.backgroundImage = ''
        win.backgroundFit = 'none'
        win.units = 'deg'
    if expInfo is not None:
        # get/measure frame rate if not already in expInfo
        if win._monitorFrameRate is None:
            win._monitorFrameRate = win.getActualFrameRate(infoMsg='Attempting to measure frame rate of screen, please wait...')
        expInfo['frameRate'] = win._monitorFrameRate
    win.hideMessage()
    if PILOTING:
        # show a visual indicator if we're in piloting mode
        if prefs.piloting['showPilotingIndicator']:
            win.showPilotingIndicator()
        # always show the mouse in piloting mode
        if prefs.piloting['forceMouseVisible']:
            win.mouseVisible = True
    
    return win


def setupDevices(expInfo, thisExp, win):
    """
    Setup whatever devices are available (mouse, keyboard, speaker, eyetracker, etc.) and add them to 
    the device manager (deviceManager)
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window in which to run this experiment.
    Returns
    ==========
    bool
        True if completed successfully.
    """
    # --- Setup input devices ---
    ioConfig = {}
    
    # Setup iohub keyboard
    ioConfig['Keyboard'] = dict(use_keymap='psychopy')
    
    # Setup iohub experiment
    ioConfig['Experiment'] = dict(filename=thisExp.dataFileName)
    
    # Start ioHub server
    ioServer = io.launchHubServer(window=win, **ioConfig)
    
    # store ioServer object in the device manager
    deviceManager.ioServer = ioServer
    
    # create a default keyboard (e.g. to check for escape)
    if deviceManager.getDevice('defaultKeyboard') is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='iohub'
        )
    # return True if completed successfully
    return True

def pauseExperiment(thisExp, win=None, timers=[], currentRoutine=None):
    """
    Pause this experiment, preventing the flow from advancing to the next routine until resumed.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    timers : list, tuple
        List of timers to reset once pausing is finished.
    currentRoutine : psychopy.data.Routine
        Current Routine we are in at time of pausing, if any. This object tells PsychoPy what Components to pause/play/dispatch.
    """
    # if we are not paused, do nothing
    if thisExp.status != PAUSED:
        return
    
    # start a timer to figure out how long we're paused for
    pauseTimer = core.Clock()
    # pause any playback components
    if currentRoutine is not None:
        for comp in currentRoutine.getPlaybackComponents():
            comp.pause()
    # make sure we have a keyboard
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        defaultKeyboard = deviceManager.addKeyboard(
            deviceClass='keyboard',
            deviceName='defaultKeyboard',
            backend='ioHub',
        )
    # run a while loop while we wait to unpause
    while thisExp.status == PAUSED:
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=['escape']):
            endExperiment(thisExp, win=win)
        # dispatch messages on response components
        if currentRoutine is not None:
            for comp in currentRoutine.getDispatchComponents():
                comp.device.dispatchMessages()
        # sleep 1ms so other threads can execute
        clock.time.sleep(0.001)
    # if stop was requested while paused, quit
    if thisExp.status == FINISHED:
        endExperiment(thisExp, win=win)
    # resume any playback components
    if currentRoutine is not None:
        for comp in currentRoutine.getPlaybackComponents():
            comp.play()
    # reset any timers
    for timer in timers:
        timer.addTime(-pauseTimer.getTime())


def run(expInfo, thisExp, win, globalClock=None, thisSession=None):
    """
    Run the experiment flow.
    
    Parameters
    ==========
    expInfo : dict
        Information about this experiment, created by the `setupExpInfo` function.
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    psychopy.visual.Window
        Window in which to run this experiment.
    globalClock : psychopy.core.clock.Clock or None
        Clock to get global time from - supply None to make a new one.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    # enter 'rush' mode (raise CPU priority)
    if not PILOTING:
        core.rush(enable=True)
    # mark experiment as started
    thisExp.status = STARTED
    # update experiment info
    expInfo['date'] = data.getDateStr()
    expInfo['expName'] = expName
    expInfo['expVersion'] = expVersion
    expInfo['psychopyVersion'] = psychopyVersion
    # make sure window is set to foreground to prevent losing focus
    win.winHandle.activate()
    # make sure variables created by exec are available globally
    exec = environmenttools.setExecEnvironment(globals())
    # get device handles from dict of input devices
    ioServer = deviceManager.ioServer
    # get/create a default keyboard (e.g. to check for escape)
    defaultKeyboard = deviceManager.getDevice('defaultKeyboard')
    if defaultKeyboard is None:
        deviceManager.addDevice(
            deviceClass='keyboard', deviceName='defaultKeyboard', backend='ioHub'
        )
    eyetracker = deviceManager.getDevice('eyetracker')
    # make sure we're running in the directory for this experiment
    os.chdir(_thisDir)
    # get filename from ExperimentHandler for convenience
    filename = thisExp.dataFileName
    frameTolerance = 0.001  # how close to onset before 'same' frame
    endExpNow = False  # flag for 'escape' or other condition => quit the exp
    # get frame duration from frame rate in expInfo
    if 'frameRate' in expInfo and expInfo['frameRate'] is not None:
        frameDur = 1.0 / round(expInfo['frameRate'])
    else:
        frameDur = 1.0 / 60.0  # could not measure, so guess
    
    # Start Code - component code to be run after the window creation
    
    # --- Initialize components for Routine "setup" ---
    # Run 'Begin Experiment' code from define_vars
    core.rush(True, realtime=True)
    
    if expInfo['EEG'] == '1':
        print('Set port')
        from psychopy import parallel
        try:
            port_address = int(expInfo['port address'], 0)   # '0' base allows both hex & decimal
            port = parallel.ParallelPort(address=port_address)
            port.setData(0)
            print(f"Parallel port initialized at {hex(port_address)}")
        except ValueError:
            print(f"Invalid port address: {expInfo['port address']}")
        except Exception as e:
            print(f"Error initializing parallel port: {e}")
    
    #trigdic = {
    #    'start': '1', 
    #    'notStart': '0', 
    #    '1+2': '1100', 
    #    '1+3': '1010', 
    #    '2+4': '0101', 
    #    '3+4': '0011', 
    #    'left': '0', 
    #    'right': '1', 
    #    'notResp': '0', 
    #    'resp': '1'
    #    }
    
    trigdic = {
        'start': '11',  
        'resp': '00', 
        'calib': '10',
        'trianing': '1',
        'experiment': '0',
        '1+2': '1100', 
        '1+3': '1010', 
        '2+4': '0101', 
        '3+4': '0011', 
        '0': '0001',
        '1': '0010',
        '2': '0011',
        '3': '0100',
        '4': '0101',
        '5': '0110',
        '6': '0111'
        }
    
    
    trigN = str()
    
    def sendTrigger(t, trigN, EEG, trigDur):
        if str(EEG) == '1':
            if t < trigDur and t > 0:  # send trigger for 50 ms and do not send the trigger before next flip time
                port.setData(int(trigN, 2)) 
            else:
                port.setData(0)
    
    instrFile = 'translations/translation_'+ expInfo['translation'].lower()+'.xlsx'
    
    # Multiplier that scales the presentation times if the refresh rate is not 60
    frameConst = int(expInfo['refreshRate'])/60 
    
    # This is to sync at 533.3(3) ms after the start of the presentation
    startFrame = 390 * frameConst 
    
    # Waiting response after each presentation for that amount of time 
    # Meaning that the "blank period" is this + whatever gets defined in iti routine
    waitRespTime = 1 
    
    # NB! Actual stimulus duration is stimDur - waitRespTime 
    stimDur = 3.1 + waitRespTime # stimDur = 4.1 # (e.g. 186 + 60 frames with 60 Hz)
    
    k = 0.8 # Just a scaler for rescaling the sizes 
    boxSize = (2.5, 3.2) # 2.5, 3.2
    symbolSize = (2.5*k, 3.2*k) # Should scale down compared to boxes (if tranparent)?
    
    # Distance from centre to the edge of the stimuli in deg
    ecc = [3, 8] # 5.25, 10.25 x (4 or 9) + (width of the stimulus (2.5 deg) / 2)
    # List of box positions 
    xys = [(ecc[1]*-1,0), (ecc[0]*-1,0),(ecc[0],0), (ecc[1],0)] 
    
    # Array of place holder boxes to be presented on screen
    rects = visual.ElementArrayStim(
        win, name = 'rects', units='deg', 
        fieldPos=(0.0, 0.0), fieldSize=(3, 4), fieldShape='square', 
        nElements=4, sizes= boxSize, xys=xys, 
        colors=([1.0, 1.0, 1.0]) , colorSpace='rgb', opacities=1, oris=0, 
        sfs=0, contrs=[1, 1,1,1], phases=0, elementTex='sqr',
        elementMask=None, texRes=48, interpolate=True, 
        autoLog=None, maskParams=None
    )
    
    # Make a keyboard object
    kb = keyboard.Keyboard()
    
    # List of images used
    imageArray = [
        'stimuli/rect_ur.png',
        'stimuli/rect_dr.png',
        'stimuli/rect_dl.png',
        'stimuli/rect_ul.png',
        'stimuli/rect_target.png'
    ]
    
    # 4x as many target stimuli in the set should give 
    # (random sampling without replacement) aprox. 70% trials with 1-3 target pairs
    randImage = [0,1,2,3,4,4,4,4] # Each number represents a symbol (4 == target)
    shuffle(randImage) # Shuffle the pool for the first round
    
    # A set with no targets (to replace targets if minimum interval of 905 ms have not eceeded after presenting the last target pair)
    nonTargetSet = imageArray[0:4] 
    
    # This helps with data entries (see the code component in the trial routine)
    def addData(rt, isPair, trials, nrOfEntries, accuracy):
    #    if nrOfEntries > 0:
    #        thisExp.nextEntry()
        thisExp.nextEntry()
        thisExp.addData('hand', save_hand)
        thisExp.addData('RT', rt) 
        thisExp.addData('targetPair', isPair)
        thisExp.addData('absNumOfTrials', trials)
        thisExp.addData('accuracy', accuracy)
    
    # Frame flipping schedule
    flipAfterOriginal = [4, 7, 3, 5] # Flip after that many frames (on 60 hz monitor)
    flipAfterEvery = [element * frameConst for element in flipAfterOriginal] # Multiply with the scaler
    
    # Chec if opacity of the box needs to be turned up (see the code component in the trial routine)
    def checkOpaStatus(opas, flipAfterEvery, frameNow):
        for opai in range(0, len(opas)):
            if frameNow % flipAfterEvery[opai] == 0:
                opas[opai] = 1
        return opas
    
    # Change the symbol in that many frames
    symShowFrames = 11*frameConst
    waitNextPairFrames = 54*frameConst # waitNextPairFrames = 54*frameConst # (aprox 900 ms, 905 in original)
    
    # Load the table using PsychoPy's built-in function
    task_data = data.importConditions(instrFile)
    #print(task_data)
    # Initialize the dictionary
    task_texts = {}
    
    # Populate the dictionary from the loaded table
    for row in task_data:
        key = row['key'].strip()
        message = row['translated_message'].strip()
        task_texts[key] = message
    
    # For specific use (e.g., hands in translation)
    hands = ['left', 'right']
    hands_in_translation = [task_texts['left_translation'], task_texts['right_translation']]
    
    # Example: Counting keys starting with 'general_intro'
    partial_key = 'general_intro'
    # Initialize count
    count = 0
    
    # Iterate through the keys and count matches
    for key in task_texts.keys():
        if key.startswith(partial_key):
            count += 1
    
    if int(expInfo['testRun']):
        presentInstructionsFor = 0.5
    else:
        presentInstructionsFor = 1e6  # Effectively "infinite" 
    
    vol = float(expInfo['volume'])
    intro_text = visual.TextStim(win=win, name='intro_text',
        text='',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=1.0, wrapWidth=30.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    intro_key = keyboard.Keyboard(deviceName='defaultKeyboard')
    
    # --- Initialize components for Routine "block_intro" ---
    block_text = visual.TextStim(win=win, name='block_text',
        text='',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=1.0, wrapWidth=30.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    block_key_resp = keyboard.Keyboard(deviceName='defaultKeyboard')
    
    # --- Initialize components for Routine "cond_setup" ---
    # Run 'Begin Experiment' code from cond_setup_code
    eog_calib = 1
    
    def make_positions(layout):
        if layout == "x":
            return [(float(x), 0.0) for x in x_points]
        elif layout == "grid":
            xs = np.linspace(-x_max, x_max, nx)
            ys = np.linspace(-y_max, y_max, ny)
            pos = [(float(x), float(y)) for y in ys[::-1] for x in xs]
            return pos
        else:
            raise ValueError("Unknown layout: " + str(layout))
    calib_dot = visual.ShapeStim(
        win=win, name='calib_dot',units='deg', 
        size=(0.5, 0.5), vertices='circle',
        ori=0.0, pos=(0, 0), draggable=False, anchor='center',
        lineWidth=2.0,
        colorSpace='rgb', lineColor='blue', fillColor='red',
        opacity=None, depth=-1.0, interpolate=True)
    
    # --- Initialize components for Routine "cond_intro" ---
    text_block_above = visual.TextStim(win=win, name='text_block_above',
        text='',
        font='Open Sans',
        pos=(0, 5), draggable=False, height=2.0, wrapWidth=30.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    key_block_intro = keyboard.Keyboard(deviceName='defaultKeyboard')
    text_block_below = visual.TextStim(win=win, name='text_block_below',
        text='',
        font='Open Sans',
        pos=(0, -5), draggable=False, height=2.0, wrapWidth=30.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-3.0);
    please_fixate = visual.TextStim(win=win, name='please_fixate',
        text='',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=0.75, wrapWidth=6.0, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-4.0);
    
    # --- Initialize components for Routine "iti" ---
    # Run 'Begin Experiment' code from iti_code
    absNumOfTrials = 0 # Keeps track of the number of trials 
    fix_iti = visual.TextStim(win=win, name='fix_iti',
        text='+',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-1.0);
    
    # --- Initialize components for Routine "trial" ---
    # Run 'Begin Experiment' code from presentStim
    import copy
    
    # initialize 'Speaker'
    deviceManager.addDevice(
        deviceName='Speaker',
        index=6.0,
        name='Speakers (Realtek(R) Audio)',
        deviceClass='psychopy.hardware.speaker.SpeakerDevice',
        resample=True,
        latencyClass=1,
    )
    mySound = sound.Sound('A', stereo=True, hamming=True, secs=0.180, speaker='Speaker', name='mySound')
    image_a = visual.ImageStim(
        win=win,
        name='image_a', units='deg', 
        image='stimuli/rect_target.png', mask=None, anchor='center',
        ori=0.0, pos=[xys[0]], draggable=False, size=symbolSize,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-1.0)
    image_b = visual.ImageStim(
        win=win,
        name='image_b', units='deg', 
        image='stimuli/rect_target.png', mask=None, anchor='center',
        ori=0.0, pos=[xys[1]], draggable=False, size=symbolSize,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-2.0)
    image_c = visual.ImageStim(
        win=win,
        name='image_c', units='deg', 
        image='stimuli/rect_target.png', mask=None, anchor='center',
        ori=0.0, pos=[xys[2]], draggable=False, size=symbolSize,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-3.0)
    image_d = visual.ImageStim(
        win=win,
        name='image_d', units='deg', 
        image='stimuli/rect_target.png', mask=None, anchor='center',
        ori=0.0, pos=[xys[3]], draggable=False, size=symbolSize,
        color=[1,1,1], colorSpace='rgb', opacity=None,
        flipHoriz=False, flipVert=False,
        texRes=128.0, interpolate=True, depth=-4.0)
    fix_trial = visual.TextStim(win=win, name='fix_trial',
        text='+',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=-5.0);
    
    # --- Initialize components for Routine "outro" ---
    outro_text = visual.TextStim(win=win, name='outro_text',
        text='',
        font='Open Sans',
        pos=(0, 0), draggable=False, height=1.0, wrapWidth=None, ori=0.0, 
        color='white', colorSpace='rgb', opacity=None, 
        languageStyle='LTR',
        depth=0.0);
    
    # create some handy timers
    
    # global clock to track the time since experiment started
    if globalClock is None:
        # create a clock if not given one
        globalClock = core.Clock()
    if isinstance(globalClock, str):
        # if given a string, make a clock accoridng to it
        if globalClock == 'float':
            # get timestamps as a simple value
            globalClock = core.Clock(format='float')
        elif globalClock == 'iso':
            # get timestamps in ISO format
            globalClock = core.Clock(format='%Y-%m-%d_%H:%M:%S.%f%z')
        else:
            # get timestamps in a custom format
            globalClock = core.Clock(format=globalClock)
    if ioServer is not None:
        ioServer.syncClock(globalClock)
    logging.setDefaultClock(globalClock)
    if eyetracker is not None:
        eyetracker.enableEventReporting()
    # routine timer to track time remaining of each (possibly non-slip) routine
    routineTimer = core.Clock()
    win.flip()  # flip window to reset last flip timer
    # store the exact time the global clock started
    expInfo['expStart'] = data.getDateStr(
        format='%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6
    )
    
    # set up handler to look after randomisation of conditions etc
    intro = data.TrialHandler2(
        name='intro',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions(
        instrFile, 
        selection=range(0,count)
    )
    , 
        seed=None, 
        isTrials=False, 
    )
    thisExp.addLoop(intro)  # add the loop to the experiment
    thisIntro = intro.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisIntro.rgb)
    if thisIntro != None:
        for paramName in thisIntro:
            globals()[paramName] = thisIntro[paramName]
    
    for thisIntro in intro:
        intro.status = STARTED
        if hasattr(thisIntro, 'status'):
            thisIntro.status = STARTED
        currentLoop = intro
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisIntro.rgb)
        if thisIntro != None:
            for paramName in thisIntro:
                globals()[paramName] = thisIntro[paramName]
        
        # --- Prepare to start Routine "setup" ---
        # create an object to store info about Routine setup
        setup = data.Routine(
            name='setup',
            components=[intro_text, intro_key],
        )
        setup.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        intro_text.setText(translated_message)
        # create starting attributes for intro_key
        intro_key.keys = []
        intro_key.rt = []
        _intro_key_allKeys = []
        # store start times for setup
        setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        setup.tStart = globalClock.getTime(format='float')
        setup.status = STARTED
        thisExp.addData('setup.started', setup.tStart)
        setup.maxDuration = None
        # keep track of which components have finished
        setupComponents = setup.components
        for thisComponent in setup.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "setup" ---
        thisExp.currentRoutine = setup
        setup.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # if trial has changed, end Routine now
            if hasattr(thisIntro, 'status') and thisIntro.status == STOPPING:
                continueRoutine = False
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *intro_text* updates
            
            # if intro_text is starting this frame...
            if intro_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                intro_text.frameNStart = frameN  # exact frame index
                intro_text.tStart = t  # local t and not account for scr refresh
                intro_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(intro_text, 'tStartRefresh')  # time at next scr refresh
                # update status
                intro_text.status = STARTED
                intro_text.setAutoDraw(True)
            
            # if intro_text is active this frame...
            if intro_text.status == STARTED:
                # update params
                pass
            
            # *intro_key* updates
            waitOnFlip = False
            
            # if intro_key is starting this frame...
            if intro_key.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                intro_key.frameNStart = frameN  # exact frame index
                intro_key.tStart = t  # local t and not account for scr refresh
                intro_key.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(intro_key, 'tStartRefresh')  # time at next scr refresh
                # update status
                intro_key.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(intro_key.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(intro_key.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if intro_key.status == STARTED and not waitOnFlip:
                theseKeys = intro_key.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                _intro_key_allKeys.extend(theseKeys)
                if len(_intro_key_allKeys):
                    intro_key.keys = _intro_key_allKeys[-1].name  # just the last key pressed
                    intro_key.rt = _intro_key_allKeys[-1].rt
                    intro_key.duration = _intro_key_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer, globalClock], 
                    currentRoutine=setup,
                )
                # skip the frame we paused on
                continue
            
            # has a Component requested the Routine to end?
            if not continueRoutine:
                setup.forceEnded = routineForceEnded = True
            # has the Routine been forcibly ended?
            if setup.forceEnded or routineForceEnded:
                break
            # has every Component finished?
            continueRoutine = False
            for thisComponent in setup.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "setup" ---
        for thisComponent in setup.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for setup
        setup.tStop = globalClock.getTime(format='float')
        setup.tStopRefresh = tThisFlipGlobal
        thisExp.addData('setup.stopped', setup.tStop)
        # Run 'End Routine' code from define_vars
        thisExp.addData('translated_message',expInfo['translation'].lower()+'.xlsx')
        # check responses
        if intro_key.keys in ['', [], None]:  # No response was made
            intro_key.keys = None
        intro.addData('intro_key.keys',intro_key.keys)
        if intro_key.keys != None:  # we had a response
            intro.addData('intro_key.rt', intro_key.rt)
            intro.addData('intro_key.duration', intro_key.duration)
        # the Routine "setup" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        # mark thisIntro as finished
        if hasattr(thisIntro, 'status'):
            thisIntro.status = FINISHED
        # if awaiting a pause, pause now
        if intro.status == PAUSED:
            thisExp.status = PAUSED
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[globalClock], 
            )
            # once done pausing, restore running status
            intro.status = STARTED
    # completed 1.0 repeats of 'intro'
    intro.status = FINISHED
    
    
    # set up handler to look after randomisation of conditions etc
    blocks = data.TrialHandler2(
        name='blocks',
        nReps=1.0, 
        method='sequential', 
        extraInfo=expInfo, 
        originPath=-1, 
        trialList=data.importConditions('blocks.xlsx'), 
        seed=None, 
        isTrials=False, 
    )
    thisExp.addLoop(blocks)  # add the loop to the experiment
    thisBlock = blocks.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock:
            globals()[paramName] = thisBlock[paramName]
    
    for thisBlock in blocks:
        blocks.status = STARTED
        if hasattr(thisBlock, 'status'):
            thisBlock.status = STARTED
        currentLoop = blocks
        thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
        # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
        if thisBlock != None:
            for paramName in thisBlock:
                globals()[paramName] = thisBlock[paramName]
        
        # --- Prepare to start Routine "block_intro" ---
        # create an object to store info about Routine block_intro
        block_intro = data.Routine(
            name='block_intro',
            components=[block_text, block_key_resp],
        )
        block_intro.status = NOT_STARTED
        continueRoutine = True
        # update component parameters for each repeat
        # Run 'Begin Routine' code from block_code
        if 'training' in condFile:
            bText = task_texts['bText_training'] # 'Following are the training trials \n\n Press "space" to begin...'
            selectRows = list(range(0,1))
        else:
            bText = task_texts['bText_experiment'] # 'Following are the experimental trials \n\n Press "space" to begin...'
            selectRows = ''
        block_text.setText(bText)
        # create starting attributes for block_key_resp
        block_key_resp.keys = []
        block_key_resp.rt = []
        _block_key_resp_allKeys = []
        # store start times for block_intro
        block_intro.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
        block_intro.tStart = globalClock.getTime(format='float')
        block_intro.status = STARTED
        block_intro.maxDuration = presentInstructionsFor
        # keep track of which components have finished
        block_introComponents = block_intro.components
        for thisComponent in block_intro.components:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "block_intro" ---
        thisExp.currentRoutine = block_intro
        block_intro.forceEnded = routineForceEnded = not continueRoutine
        while continueRoutine:
            # if trial has changed, end Routine now
            if hasattr(thisBlock, 'status') and thisBlock.status == STOPPING:
                continueRoutine = False
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # is it time to end the Routine? (based on local clock)
            if tThisFlip > block_intro.maxDuration-frameTolerance:
                block_intro.maxDurationReached = True
                continueRoutine = False
            
            # *block_text* updates
            
            # if block_text is starting this frame...
            if block_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_text.frameNStart = frameN  # exact frame index
                block_text.tStart = t  # local t and not account for scr refresh
                block_text.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_text, 'tStartRefresh')  # time at next scr refresh
                # update status
                block_text.status = STARTED
                block_text.setAutoDraw(True)
            
            # if block_text is active this frame...
            if block_text.status == STARTED:
                # update params
                pass
            
            # *block_key_resp* updates
            waitOnFlip = False
            
            # if block_key_resp is starting this frame...
            if block_key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                block_key_resp.frameNStart = frameN  # exact frame index
                block_key_resp.tStart = t  # local t and not account for scr refresh
                block_key_resp.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(block_key_resp, 'tStartRefresh')  # time at next scr refresh
                # update status
                block_key_resp.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(block_key_resp.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(block_key_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if block_key_resp.status == STARTED and not waitOnFlip:
                theseKeys = block_key_resp.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                _block_key_resp_allKeys.extend(theseKeys)
                if len(_block_key_resp_allKeys):
                    block_key_resp.keys = _block_key_resp_allKeys[-1].name  # just the last key pressed
                    block_key_resp.rt = _block_key_resp_allKeys[-1].rt
                    block_key_resp.duration = _block_key_resp_allKeys[-1].duration
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if defaultKeyboard.getKeys(keyList=["escape"]):
                thisExp.status = FINISHED
            if thisExp.status == FINISHED or endExpNow:
                endExperiment(thisExp, win=win)
                return
            # pause experiment here if requested
            if thisExp.status == PAUSED:
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[routineTimer, globalClock], 
                    currentRoutine=block_intro,
                )
                # skip the frame we paused on
                continue
            
            # has a Component requested the Routine to end?
            if not continueRoutine:
                block_intro.forceEnded = routineForceEnded = True
            # has the Routine been forcibly ended?
            if block_intro.forceEnded or routineForceEnded:
                break
            # has every Component finished?
            continueRoutine = False
            for thisComponent in block_intro.components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "block_intro" ---
        for thisComponent in block_intro.components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # store stop times for block_intro
        block_intro.tStop = globalClock.getTime(format='float')
        block_intro.tStopRefresh = tThisFlipGlobal
        # check responses
        if block_key_resp.keys in ['', [], None]:  # No response was made
            block_key_resp.keys = None
        blocks.addData('block_key_resp.keys',block_key_resp.keys)
        if block_key_resp.keys != None:  # we had a response
            blocks.addData('block_key_resp.rt', block_key_resp.rt)
            blocks.addData('block_key_resp.duration', block_key_resp.duration)
        # the Routine "block_intro" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        conditions = data.TrialHandler2(
            name='conditions',
            nReps=1.0, 
            method='random', 
            extraInfo=expInfo, 
            originPath=-1, 
            trialList=data.importConditions(condFile), 
            seed=None, 
            isTrials=True, 
        )
        thisExp.addLoop(conditions)  # add the loop to the experiment
        thisCondition = conditions.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisCondition.rgb)
        if thisCondition != None:
            for paramName in thisCondition:
                globals()[paramName] = thisCondition[paramName]
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        
        for thisCondition in conditions:
            conditions.status = STARTED
            if hasattr(thisCondition, 'status'):
                thisCondition.status = STARTED
            currentLoop = conditions
            thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            # abbreviate parameter names if possible (e.g. rgb = thisCondition.rgb)
            if thisCondition != None:
                for paramName in thisCondition:
                    globals()[paramName] = thisCondition[paramName]
            
            # --- Prepare to start Routine "cond_setup" ---
            # create an object to store info about Routine cond_setup
            cond_setup = data.Routine(
                name='cond_setup',
                components=[calib_dot],
            )
            cond_setup.status = NOT_STARTED
            continueRoutine = True
            # update component parameters for each repeat
            # Run 'Begin Routine' code from cond_setup_code
            if isTraining:
                nTrials = 10 # trials per each row in the condition table
            else:
                nTrials = 40
            
            if eog_calib and conditions.thisN % 4 == 0:
                # Choose calibration layout
                calib_layout = "x"     # "x" (horizontal only) or "grid" (2D)
            
                x_points = [-ecc[1], -ecc[0], -1, 0, 1, ecc[0], ecc[1]]   # or build from ecc + extras
            
                # If you ever need a 2D grid later
                nx, ny = 5, 5
                x_max, y_max = max(x_points), max(x_points)   # deg (set once)
            
                calibration_positions = make_positions(calib_layout)
                shuffle(calibration_positions)
                # ensure center first (optional)
                calibration_positions = [p for p in calibration_positions if p != (0.0, 0.0)]
                calibration_positions.insert(0, (0.0, 0.0))
            
                # timing/state as before...
                fixation_hold_time = 0.5
                sample_time_default = 1.5
                sample_time_first = 2.0
            
                current_calib_index = 0
                phase = "hold"
                phase_start = 0.0
                calibration_done = False
            
                calib_dot.opacity = 1.0
                calib_dot.pos = calibration_positions[current_calib_index]
            else:
                continueRoutine = False
            
            
            
            # store start times for cond_setup
            cond_setup.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
            cond_setup.tStart = globalClock.getTime(format='float')
            cond_setup.status = STARTED
            cond_setup.maxDuration = None
            # keep track of which components have finished
            cond_setupComponents = cond_setup.components
            for thisComponent in cond_setup.components:
                thisComponent.tStart = None
                thisComponent.tStop = None
                thisComponent.tStartRefresh = None
                thisComponent.tStopRefresh = None
                if hasattr(thisComponent, 'status'):
                    thisComponent.status = NOT_STARTED
            # reset timers
            t = 0
            _timeToFirstFrame = win.getFutureFlipTime(clock="now")
            frameN = -1
            
            # --- Run Routine "cond_setup" ---
            thisExp.currentRoutine = cond_setup
            cond_setup.forceEnded = routineForceEnded = not continueRoutine
            while continueRoutine:
                # if trial has changed, end Routine now
                if hasattr(thisCondition, 'status') and thisCondition.status == STOPPING:
                    continueRoutine = False
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from cond_setup_code
                current_point = x_points.index(calibration_positions[current_calib_index][0])
                if t > 2:
                    trigN = '1'+ str(isTraining) + trigdic['calib'] + trigdic[str(current_point)]
                    #print(x_points.index(calibration_positions[current_calib_index][0]))
                    # --- Key handling ---
                    keys = event.getKeys()
                    if 'space' in keys and eog_calib:
                        # restart calibration
                        current_calib_index = 0
                        phase = "hold"
                        phase_start = t
                        calibration_done = False
                        calib_dot.opacity = 1.0
                        calib_dot.pos = calibration_positions[current_calib_index]
                    # --- Calibration logic ---
                    if eog_calib and not calibration_done:
                        # durations
                        sample_duration = sample_time_first if current_calib_index == 0 else sample_time_default
                
                        # always keep dot at the current target
                        calib_dot.pos = calibration_positions[current_calib_index]
                
                        elapsed = t - phase_start
                
                        if phase == "hold":
                            # steady dot
                            calib_dot.lineWidth = 2.0
                
                            if elapsed >= fixation_hold_time:
                                phase = "sample"
                                phase_start = t
                
                        elif phase == "sample":
                            # pulse lineWidth a bit for visibility (purely cosmetic)
                            # smooth 0..1 progress through the sample window:
                            p = min(1.0, elapsed / sample_duration)
                            calib_dot.lineWidth = 5.0 + 5.0 * np.sin(2*np.pi * p -np.pi/4) 
                
                            if elapsed >= sample_duration:
                                # advance to next point
                                current_calib_index += 1
                                phase = "hold"
                                phase_start = t
                
                                if current_calib_index >= len(calibration_positions):
                                    calibration_done = True
                                    calib_dot.opacity = 0.0  # hide dot when done
                                    continueRoutine = False
                
                        calib_dot.draw()
                
                # *calib_dot* updates
                
                # if calib_dot is starting this frame...
                if calib_dot.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    calib_dot.frameNStart = frameN  # exact frame index
                    calib_dot.tStart = t  # local t and not account for scr refresh
                    calib_dot.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(calib_dot, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'calib_dot.started')
                    # update status
                    calib_dot.status = STARTED
                    calib_dot.setAutoDraw(True)
                
                # if calib_dot is active this frame...
                if calib_dot.status == STARTED:
                    # update params
                    pass
                
                # check for quit (typically the Esc key)
                if defaultKeyboard.getKeys(keyList=["escape"]):
                    thisExp.status = FINISHED
                if thisExp.status == FINISHED or endExpNow:
                    endExperiment(thisExp, win=win)
                    return
                # pause experiment here if requested
                if thisExp.status == PAUSED:
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[routineTimer, globalClock], 
                        currentRoutine=cond_setup,
                    )
                    # skip the frame we paused on
                    continue
                
                # has a Component requested the Routine to end?
                if not continueRoutine:
                    cond_setup.forceEnded = routineForceEnded = True
                # has the Routine been forcibly ended?
                if cond_setup.forceEnded or routineForceEnded:
                    break
                # has every Component finished?
                continueRoutine = False
                for thisComponent in cond_setup.components:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "cond_setup" ---
            for thisComponent in cond_setup.components:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # store stop times for cond_setup
            cond_setup.tStop = globalClock.getTime(format='float')
            cond_setup.tStopRefresh = tThisFlipGlobal
            # the Routine "cond_setup" was not non-slip safe, so reset the non-slip timer
            routineTimer.reset()
            
            # set up handler to look after randomisation of conditions etc
            trials = data.TrialHandler2(
                name='trials',
                nReps=nTrials, 
                method='random', 
                extraInfo=expInfo, 
                originPath=-1, 
                trialList=[None], 
                seed=None, 
                isTrials=True, 
            )
            thisExp.addLoop(trials)  # add the loop to the experiment
            thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
            # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
            if thisTrial != None:
                for paramName in thisTrial:
                    globals()[paramName] = thisTrial[paramName]
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            
            for thisTrial in trials:
                trials.status = STARTED
                if hasattr(thisTrial, 'status'):
                    thisTrial.status = STARTED
                currentLoop = trials
                thisExp.timestampOnFlip(win, 'thisRow.t', format=globalClock.format)
                if thisSession is not None:
                    # if running in a Session with a Liaison client, send data up to now
                    thisSession.sendExperimentData()
                # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
                if thisTrial != None:
                    for paramName in thisTrial:
                        globals()[paramName] = thisTrial[paramName]
                
                # --- Prepare to start Routine "cond_intro" ---
                # create an object to store info about Routine cond_intro
                cond_intro = data.Routine(
                    name='cond_intro',
                    components=[text_block_above, key_block_intro, text_block_below, please_fixate],
                )
                cond_intro.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # Run 'Begin Routine' code from hands
                halfOftrials = int(nTrials/2)
                
                # Randomize hands before each block
                if not trials.thisN:
                    hand_list = [0,1]
                    shuffle(hand_list)
                
                if trials.thisN > halfOftrials-1:
                    current_hand = hands_in_translation[hand_list[0]] # for presenting
                    save_hand = hands[hand_list[0]] # for saving
                else:
                    current_hand = hands_in_translation[hand_list[1]]
                    save_hand = hands[hand_list[1]]
                
                # Define numbers to draw
                numbers = ['1', '2', '3', '4']
                
                # Create TextStim objects for each number
                text_stims = []
                for num, pos in zip(numbers, xys):
                    text_stim = visual.TextStim(
                        win=win,
                        text=num,
                        pos=pos,
                        color='black',
                        height=1.0  # Adjust as needed
                    )
                    text_stims.append(text_stim)
                text_block_above.setText(task_texts['cond_text_attend'] + ' ' + attendCond)
                # create starting attributes for key_block_intro
                key_block_intro.keys = []
                key_block_intro.rt = []
                _key_block_intro_allKeys = []
                text_block_below.setText(task_texts['cond_text_hand'] + ' ' + current_hand)
                please_fixate.setText(task_texts['gaze'])
                # store start times for cond_intro
                cond_intro.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                cond_intro.tStart = globalClock.getTime(format='float')
                cond_intro.status = STARTED
                cond_intro.maxDuration = presentInstructionsFor
                # skip Routine cond_intro if its 'Skip if' condition is True
                cond_intro.skipped = continueRoutine and not (trials.thisN != int(nTrials/2) and trials.thisN != 0)
                continueRoutine = cond_intro.skipped
                # keep track of which components have finished
                cond_introComponents = cond_intro.components
                for thisComponent in cond_intro.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "cond_intro" ---
                thisExp.currentRoutine = cond_intro
                cond_intro.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine:
                    # if trial has changed, end Routine now
                    if hasattr(thisTrial, 'status') and thisTrial.status == STOPPING:
                        continueRoutine = False
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    # is it time to end the Routine? (based on local clock)
                    if tThisFlip > cond_intro.maxDuration-frameTolerance:
                        cond_intro.maxDurationReached = True
                        continueRoutine = False
                    # Run 'Each Frame' code from hands
                    rects.draw() # Present the boxes
                    
                    # Draw all numbers
                    for stim in text_stims:
                        stim.draw()
                    
                    # *text_block_above* updates
                    
                    # if text_block_above is starting this frame...
                    if text_block_above.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        text_block_above.frameNStart = frameN  # exact frame index
                        text_block_above.tStart = t  # local t and not account for scr refresh
                        text_block_above.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(text_block_above, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        text_block_above.status = STARTED
                        text_block_above.setAutoDraw(True)
                    
                    # if text_block_above is active this frame...
                    if text_block_above.status == STARTED:
                        # update params
                        pass
                    
                    # *key_block_intro* updates
                    waitOnFlip = False
                    
                    # if key_block_intro is starting this frame...
                    if key_block_intro.status == NOT_STARTED and tThisFlip >= 0.5-frameTolerance:
                        # keep track of start time/frame for later
                        key_block_intro.frameNStart = frameN  # exact frame index
                        key_block_intro.tStart = t  # local t and not account for scr refresh
                        key_block_intro.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(key_block_intro, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        key_block_intro.status = STARTED
                        # keyboard checking is just starting
                        waitOnFlip = True
                        win.callOnFlip(key_block_intro.clock.reset)  # t=0 on next screen flip
                        win.callOnFlip(key_block_intro.clearEvents, eventType='keyboard')  # clear events on next screen flip
                    if key_block_intro.status == STARTED and not waitOnFlip:
                        theseKeys = key_block_intro.getKeys(keyList=['space'], ignoreKeys=["escape"], waitRelease=False)
                        _key_block_intro_allKeys.extend(theseKeys)
                        if len(_key_block_intro_allKeys):
                            key_block_intro.keys = _key_block_intro_allKeys[-1].name  # just the last key pressed
                            key_block_intro.rt = _key_block_intro_allKeys[-1].rt
                            key_block_intro.duration = _key_block_intro_allKeys[-1].duration
                            # a response ends the routine
                            continueRoutine = False
                    
                    # *text_block_below* updates
                    
                    # if text_block_below is starting this frame...
                    if text_block_below.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        text_block_below.frameNStart = frameN  # exact frame index
                        text_block_below.tStart = t  # local t and not account for scr refresh
                        text_block_below.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(text_block_below, 'tStartRefresh')  # time at next scr refresh
                        # update status
                        text_block_below.status = STARTED
                        text_block_below.setAutoDraw(True)
                    
                    # if text_block_below is active this frame...
                    if text_block_below.status == STARTED:
                        # update params
                        pass
                    
                    # *please_fixate* updates
                    
                    # if please_fixate is starting this frame...
                    if please_fixate.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        please_fixate.frameNStart = frameN  # exact frame index
                        please_fixate.tStart = t  # local t and not account for scr refresh
                        please_fixate.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(please_fixate, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'please_fixate.started')
                        # update status
                        please_fixate.status = STARTED
                        please_fixate.setAutoDraw(True)
                    
                    # if please_fixate is active this frame...
                    if please_fixate.status == STARTED:
                        # update params
                        pass
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer, globalClock], 
                            currentRoutine=cond_intro,
                        )
                        # skip the frame we paused on
                        continue
                    
                    # has a Component requested the Routine to end?
                    if not continueRoutine:
                        cond_intro.forceEnded = routineForceEnded = True
                    # has the Routine been forcibly ended?
                    if cond_intro.forceEnded or routineForceEnded:
                        break
                    # has every Component finished?
                    continueRoutine = False
                    for thisComponent in cond_intro.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "cond_intro" ---
                for thisComponent in cond_intro.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for cond_intro
                cond_intro.tStop = globalClock.getTime(format='float')
                cond_intro.tStopRefresh = tThisFlipGlobal
                # check responses
                if key_block_intro.keys in ['', [], None]:  # No response was made
                    key_block_intro.keys = None
                trials.addData('key_block_intro.keys',key_block_intro.keys)
                if key_block_intro.keys != None:  # we had a response
                    trials.addData('key_block_intro.rt', key_block_intro.rt)
                    trials.addData('key_block_intro.duration', key_block_intro.duration)
                # the Routine "cond_intro" was not non-slip safe, so reset the non-slip timer
                routineTimer.reset()
                
                # --- Prepare to start Routine "iti" ---
                # create an object to store info about Routine iti
                iti = data.Routine(
                    name='iti',
                    components=[fix_iti],
                )
                iti.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # Run 'Begin Routine' code from iti_code
                shuffle(randImage)
                rects.opacities = 1 # Set the opacity back up for all the boxes 
                absNumOfTrials += 1 # Increase the trial counter by 1
                thisExp.addData('absNumOfTrials', absNumOfTrials) # Send info about trial number to the data file
                # store start times for iti
                iti.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                iti.tStart = globalClock.getTime(format='float')
                iti.status = STARTED
                thisExp.addData('iti.started', iti.tStart)
                iti.maxDuration = None
                # keep track of which components have finished
                itiComponents = iti.components
                for thisComponent in iti.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "iti" ---
                thisExp.currentRoutine = iti
                iti.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine and routineTimer.getTime() < 2.0:
                    # if trial has changed, end Routine now
                    if hasattr(thisTrial, 'status') and thisTrial.status == STOPPING:
                        continueRoutine = False
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    # Run 'Each Frame' code from iti_code
                    rects.draw() # Presents place holders (white boxes)
                    
                    # *fix_iti* updates
                    
                    # if fix_iti is starting this frame...
                    if fix_iti.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        fix_iti.frameNStart = frameN  # exact frame index
                        fix_iti.tStart = t  # local t and not account for scr refresh
                        fix_iti.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(fix_iti, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'fix_iti.started')
                        # update status
                        fix_iti.status = STARTED
                        fix_iti.setAutoDraw(True)
                    
                    # if fix_iti is active this frame...
                    if fix_iti.status == STARTED:
                        # update params
                        pass
                    
                    # if fix_iti is stopping this frame...
                    if fix_iti.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > fix_iti.tStartRefresh + 2-frameTolerance:
                            # keep track of stop time/frame for later
                            fix_iti.tStop = t  # not accounting for scr refresh
                            fix_iti.tStopRefresh = tThisFlipGlobal  # on global time
                            fix_iti.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'fix_iti.stopped')
                            # update status
                            fix_iti.status = FINISHED
                            fix_iti.setAutoDraw(False)
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer, globalClock], 
                            currentRoutine=iti,
                        )
                        # skip the frame we paused on
                        continue
                    
                    # has a Component requested the Routine to end?
                    if not continueRoutine:
                        iti.forceEnded = routineForceEnded = True
                    # has the Routine been forcibly ended?
                    if iti.forceEnded or routineForceEnded:
                        break
                    # has every Component finished?
                    continueRoutine = False
                    for thisComponent in iti.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "iti" ---
                for thisComponent in iti.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for iti
                iti.tStop = globalClock.getTime(format='float')
                iti.tStopRefresh = tThisFlipGlobal
                thisExp.addData('iti.stopped', iti.tStop)
                # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
                if iti.maxDurationReached:
                    routineTimer.addTime(-iti.maxDuration)
                elif iti.forceEnded:
                    routineTimer.reset()
                else:
                    routineTimer.addTime(-2.000000)
                
                # --- Prepare to start Routine "trial" ---
                # create an object to store info about Routine trial
                trial = data.Routine(
                    name='trial',
                    components=[image_a, image_b, image_c, image_d, fix_trial],
                )
                trial.status = NOT_STARTED
                continueRoutine = True
                # update component parameters for each repeat
                # Run 'Begin Routine' code from presentStim
                frameCount = 0
                responseGiven = False
                targetPair = False
                nrOfEntries = 0
                cond = [pos2attLeft-1, pos2attRight-1]
                counter = 0
                opacities = [0,0,0,0]
                switchTime = 0
                t2 = 0
                
                # Response window in seconds
                MIN_RT = 0.250
                MAX_RT = 1.000
                INTERMEDIATE_IDX = 3  # 4th position
                
                eventType = None         # 'hit', 'fa', or None
                eventOnset = None
                eventResponded = False
                
                # List of components
                imList = [image_a, image_b, image_c, image_d]
                
                #trigN = '1'+ trigdic['start'] + trigdic[attendCond] + trigdic[save_hand] + trigdic['notResp']
                trigN = '1'+ str(isTraining) + trigdic['start'] + trigdic[attendCond]
                
                # Make all the symbols visible
                for im in imList:
                    im.opacity = 1
                keys = kb.getKeys()
                
                allPresentedSymbols = list()
                
                if int(expInfo['testRun']):
                    RTs = [(i + 1) + random() for i in range(randint(1, 3))]
                    nDummyResponses = 0
                # store start times for trial
                trial.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
                trial.tStart = globalClock.getTime(format='float')
                trial.status = STARTED
                thisExp.addData('trial.started', trial.tStart)
                trial.maxDuration = None
                # keep track of which components have finished
                trialComponents = trial.components
                for thisComponent in trial.components:
                    thisComponent.tStart = None
                    thisComponent.tStop = None
                    thisComponent.tStartRefresh = None
                    thisComponent.tStopRefresh = None
                    if hasattr(thisComponent, 'status'):
                        thisComponent.status = NOT_STARTED
                # reset timers
                t = 0
                _timeToFirstFrame = win.getFutureFlipTime(clock="now")
                frameN = -1
                
                # --- Run Routine "trial" ---
                thisExp.currentRoutine = trial
                trial.forceEnded = routineForceEnded = not continueRoutine
                while continueRoutine:
                    # if trial has changed, end Routine now
                    if hasattr(thisTrial, 'status') and thisTrial.status == STOPPING:
                        continueRoutine = False
                    # get current time
                    t = routineTimer.getTime()
                    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                    # update/draw components on each frame
                    # Run 'Each Frame' code from presentStim
                    if not isTraining:
                        sendTrigger(t-t2, trigN, expInfo['EEG'], 0.05)
                    
                    rects.draw()  # Present the boxes
                    
                    # List of default opacities, this will be written over for one frame by checkOpaStatus if necessary
                    # This is just to keep the duration of "on frame" ~ 1000/60 ms (1000/60.8 in original)
                    if sum(opacities) and frameConst-counter <= 1:  # frameConst == 1 on 60 Hz
                        opacities = [0, 0, 0, 0]
                        counter = 0
                    else:
                        counter += 1
                    
                    # Sync at 533.3(3) ms after the start
                    frameNow = startFrame + frameN
                    
                    # Change the opacity based cycle duration in frames
                    opacities = checkOpaStatus(opacities, flipAfterEvery, frameNow)  # fAE == [4, 7, 3, 5] by default
                    
                    # Change the opacity of the boxes based on checkOpaStatus()
                    if t <= stimDur - waitRespTime:
                        rects.opacities = opacities
                    else:  # If the presentation is over set back to 1
                        rects.opacities = 1
                    
                    # -----------------------------
                    # Change symbol after every 183.3(3) ms
                    # -----------------------------
                    if frameN % symShowFrames == 0 and not t > stimDur - waitRespTime:  # symShowFrames == 11 frames (on 60 Hz)
                    
                        # Take random images from the set
                        randomSet = randImage[0:4]
                        imCount = 0
                        for im in imList:
                            im.image = imageArray[randomSet[imCount]]  # randImage == [0,1,2,3,4,4,4,4]
                            imCount += 1
                    
                        # --- ADDED: compute current on-screen state (after assignment) ---
                        is_target = [('target' in im.image) for im in imList]          ### ADDED
                        att0, att1 = cond[0], cond[1]                                  ### ADDED
                        hit_now = is_target[att0] and is_target[att1]                  ### ADDED
                    
                        # If the target is in both of the target locations then
                        if hit_now:
                            # If wait time has exceeded -> start a NEW valid target-pair epoch
                            if frameN - frameCount >= waitNextPairFrames + symShowFrames:
                    
                                # --- miss-by-timeout is handled below; do NOT add miss here (keeps new rule) ---
                    
                                eventOnset = t                                         ### CHANGED (new epoch marker)
                                eventResponded = False                                 ### CHANGED
                                frameCount = frameN                                    # (same role as before)
                    
                            else:
                                # If wait time is not over yet but random sampling gave another target pair
                                # then just change one of them randomly to non-target (OLD behavior)
                                selectRandLoc = cond[randint(0, 1)]
                                replaceRandSymb = randint(0, 4)
                                imList[selectRandLoc].image = nonTargetSet[replaceRandSymb]
                                randomSet[selectRandLoc] = replaceRandSymb             # keep logging consistent
                    
                                # --- ADDED: recompute state after breaking the pair ---
                                is_target = [('target' in im.image) for im in imList]  ### ADDED
                                hit_now = is_target[att0] and is_target[att1]          ### ADDED
                    
                        if not t > (stimDur - waitRespTime):
                            allPresentedSymbols.append(copy.deepcopy(randomSet))
                            if int(expInfo['testRun']):
                                if 'target' in imList[cond[0]].image and 'target' in imList[cond[1]].image:
                                    mySound = sound.Sound('A', octave=2, hamming=True, speaker='Speaker', secs=0.180, volume=vol)
                                    mySound.play()
                    
                        shuffle(randImage)  # Shuffle for the next round
                    
                    
                    # -----------------------------
                    # ADDED: miss-by-timeout (250–1000 ms window)
                    # -----------------------------
                    if (eventOnset is not None) and (not eventResponded) and ((t - eventOnset) > MAX_RT):   ### ADDED
                        addData('miss', True, absNumOfTrials, nrOfEntries, 0)                               ### ADDED
                        nrOfEntries += 1                                                                    ### ADDED
                        eventResponded = True                                                               ### ADDED
                        eventOnset = None                                                                   ### ADDED
                    
                    
                    # -----------------------------
                    # Check keys (keep original testRun simulation structure)
                    # -----------------------------
                    if int(expInfo['testRun']):
                        if len(RTs) > nDummyResponses and t >= RTs[nDummyResponses]:
                            simulated_key = keyboard.KeyPress(name='space', code='space', tDown=t)  # Simulate keypress
                            keys = [simulated_key]  # Mimic PsychoPy key structure
                            nDummyResponses += 1
                        else:
                            keys = 0
                    else:
                        keys = kb.getKeys()
                    
                    if keys:
                        if 'space' in keys[-1].name:
                            trigN = '1' + str(isTraining) + trigdic['resp'] + trigdic[attendCond]            ### CHANGED (your newer trig)
                            t2 = t
                            if not isTraining:
                                sendTrigger(t-t2, trigN, expInfo['EEG'], 0.05)
                    
                            # --- CHANGED: RT window relative to last valid target-pair onset ---
                            rt = None                                                                       ### CHANGED
                            if eventOnset is not None:                                                      ### CHANGED
                                rt = t - eventOnset                                                         ### CHANGED
                            in_window = (rt is not None) and (MIN_RT <= rt <= MAX_RT)                       ### CHANGED
                    
                            # --- ADDED: FA definition based on CURRENT display state (no FA epochs) ---
                            is_target = [('target' in im.image) for im in imList]                            ### ADDED
                            att0, att1 = cond[0], cond[1]                                                    ### ADDED
                            single_attended = (is_target[att0] ^ is_target[att1])                            ### ADDED
                            intermediate_target = is_target[INTERMEDIATE_IDX]                                ### ADDED
                    
                            if attendCond in ('1+2', '3+4'):                                                 ### ADDED
                                fa_now = single_attended                                                     ### ADDED
                            else:
                                fa_now = single_attended or intermediate_target                              ### ADDED
                    
                            # --- scoring ---
                            if in_window and (not eventResponded):                                           ### CHANGED
                                # HIT
                                addData(rt, True, absNumOfTrials, nrOfEntries, 1)                            ### CHANGED
                                nrOfEntries += 1                                                             ### CHANGED
                                eventResponded = True                                                        ### CHANGED
                                eventOnset = None                                                            ### CHANGED
                    
                                if isTraining:
                                    mySound = sound.Sound('A', octave=4, hamming=True, speaker='Speaker', secs=0.180, volume=vol)
                                    mySound.play()
                    
                            elif (eventOnset is None) and fa_now:                                            ### ADDED
                                # FALSE ALARM (only when NOT in an active target-pair epoch)
                                addData('false alarm', False, absNumOfTrials, nrOfEntries, 0)                ### ADDED
                                nrOfEntries += 1                                                             ### ADDED
                    
                            # else: ignore (press too early/late relative to target-pair epoch, or not an FA by definition)
                    
                        elif 'escape' in keys[-1].name:
                            core.quit()
                    
                    # Do not show the symbols if the presentation is over (response will be still recorded)
                    if t > stimDur - waitRespTime and imList[0].opacity:
                        for im in imList:
                            im.opacity = 0
                    
                    
                    # *image_a* updates
                    
                    # if image_a is starting this frame...
                    if image_a.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        image_a.frameNStart = frameN  # exact frame index
                        image_a.tStart = t  # local t and not account for scr refresh
                        image_a.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(image_a, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'image_a.started')
                        # update status
                        image_a.status = STARTED
                        image_a.setAutoDraw(True)
                    
                    # if image_a is active this frame...
                    if image_a.status == STARTED:
                        # update params
                        pass
                    
                    # if image_a is stopping this frame...
                    if image_a.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > image_a.tStartRefresh + stimDur-frameTolerance:
                            # keep track of stop time/frame for later
                            image_a.tStop = t  # not accounting for scr refresh
                            image_a.tStopRefresh = tThisFlipGlobal  # on global time
                            image_a.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'image_a.stopped')
                            # update status
                            image_a.status = FINISHED
                            image_a.setAutoDraw(False)
                    
                    # *image_b* updates
                    
                    # if image_b is starting this frame...
                    if image_b.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        image_b.frameNStart = frameN  # exact frame index
                        image_b.tStart = t  # local t and not account for scr refresh
                        image_b.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(image_b, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'image_b.started')
                        # update status
                        image_b.status = STARTED
                        image_b.setAutoDraw(True)
                    
                    # if image_b is active this frame...
                    if image_b.status == STARTED:
                        # update params
                        pass
                    
                    # if image_b is stopping this frame...
                    if image_b.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > image_b.tStartRefresh + stimDur-frameTolerance:
                            # keep track of stop time/frame for later
                            image_b.tStop = t  # not accounting for scr refresh
                            image_b.tStopRefresh = tThisFlipGlobal  # on global time
                            image_b.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'image_b.stopped')
                            # update status
                            image_b.status = FINISHED
                            image_b.setAutoDraw(False)
                    
                    # *image_c* updates
                    
                    # if image_c is starting this frame...
                    if image_c.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        image_c.frameNStart = frameN  # exact frame index
                        image_c.tStart = t  # local t and not account for scr refresh
                        image_c.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(image_c, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'image_c.started')
                        # update status
                        image_c.status = STARTED
                        image_c.setAutoDraw(True)
                    
                    # if image_c is active this frame...
                    if image_c.status == STARTED:
                        # update params
                        pass
                    
                    # if image_c is stopping this frame...
                    if image_c.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > image_c.tStartRefresh + stimDur-frameTolerance:
                            # keep track of stop time/frame for later
                            image_c.tStop = t  # not accounting for scr refresh
                            image_c.tStopRefresh = tThisFlipGlobal  # on global time
                            image_c.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'image_c.stopped')
                            # update status
                            image_c.status = FINISHED
                            image_c.setAutoDraw(False)
                    
                    # *image_d* updates
                    
                    # if image_d is starting this frame...
                    if image_d.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        image_d.frameNStart = frameN  # exact frame index
                        image_d.tStart = t  # local t and not account for scr refresh
                        image_d.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(image_d, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'image_d.started')
                        # update status
                        image_d.status = STARTED
                        image_d.setAutoDraw(True)
                    
                    # if image_d is active this frame...
                    if image_d.status == STARTED:
                        # update params
                        pass
                    
                    # if image_d is stopping this frame...
                    if image_d.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > image_d.tStartRefresh + stimDur-frameTolerance:
                            # keep track of stop time/frame for later
                            image_d.tStop = t  # not accounting for scr refresh
                            image_d.tStopRefresh = tThisFlipGlobal  # on global time
                            image_d.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'image_d.stopped')
                            # update status
                            image_d.status = FINISHED
                            image_d.setAutoDraw(False)
                    
                    # *fix_trial* updates
                    
                    # if fix_trial is starting this frame...
                    if fix_trial.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                        # keep track of start time/frame for later
                        fix_trial.frameNStart = frameN  # exact frame index
                        fix_trial.tStart = t  # local t and not account for scr refresh
                        fix_trial.tStartRefresh = tThisFlipGlobal  # on global time
                        win.timeOnFlip(fix_trial, 'tStartRefresh')  # time at next scr refresh
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'fix_trial.started')
                        # update status
                        fix_trial.status = STARTED
                        fix_trial.setAutoDraw(True)
                    
                    # if fix_trial is active this frame...
                    if fix_trial.status == STARTED:
                        # update params
                        pass
                    
                    # if fix_trial is stopping this frame...
                    if fix_trial.status == STARTED:
                        # is it time to stop? (based on global clock, using actual start)
                        if tThisFlipGlobal > fix_trial.tStartRefresh + stimDur-frameTolerance:
                            # keep track of stop time/frame for later
                            fix_trial.tStop = t  # not accounting for scr refresh
                            fix_trial.tStopRefresh = tThisFlipGlobal  # on global time
                            fix_trial.frameNStop = frameN  # exact frame index
                            # add timestamp to datafile
                            thisExp.timestampOnFlip(win, 'fix_trial.stopped')
                            # update status
                            fix_trial.status = FINISHED
                            fix_trial.setAutoDraw(False)
                    
                    # check for quit (typically the Esc key)
                    if defaultKeyboard.getKeys(keyList=["escape"]):
                        thisExp.status = FINISHED
                    if thisExp.status == FINISHED or endExpNow:
                        endExperiment(thisExp, win=win)
                        return
                    # pause experiment here if requested
                    if thisExp.status == PAUSED:
                        pauseExperiment(
                            thisExp=thisExp, 
                            win=win, 
                            timers=[routineTimer, globalClock], 
                            currentRoutine=trial,
                        )
                        # skip the frame we paused on
                        continue
                    
                    # has a Component requested the Routine to end?
                    if not continueRoutine:
                        trial.forceEnded = routineForceEnded = True
                    # has the Routine been forcibly ended?
                    if trial.forceEnded or routineForceEnded:
                        break
                    # has every Component finished?
                    continueRoutine = False
                    for thisComponent in trial.components:
                        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                            continueRoutine = True
                            break  # at least one component has not yet finished
                    
                    # refresh the screen
                    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                        win.flip()
                
                # --- Ending Routine "trial" ---
                for thisComponent in trial.components:
                    if hasattr(thisComponent, "setAutoDraw"):
                        thisComponent.setAutoDraw(False)
                # store stop times for trial
                trial.tStop = globalClock.getTime(format='float')
                trial.tStopRefresh = tThisFlipGlobal
                thisExp.addData('trial.stopped', trial.tStop)
                # Run 'End Routine' code from presentStim
                thisExp.addData('allStim', [allPresentedSymbols])
                
                # Final check at end of trial to log missed responses
                if not responseGiven:
                    if targetPair:
                        addData('miss', targetPair, absNumOfTrials, nrOfEntries, 0)  # Missed target
                    elif nrOfEntries == 0:
                        addData('correct rejection', targetPair, absNumOfTrials, nrOfEntries, 'CR')  # No target, no response
                    nrOfEntries += 1
                    responseGiven = True  # Mark as handled
                # the Routine "trial" was not non-slip safe, so reset the non-slip timer
                routineTimer.reset()
                # mark thisTrial as finished
                if hasattr(thisTrial, 'status'):
                    thisTrial.status = FINISHED
                # if awaiting a pause, pause now
                if trials.status == PAUSED:
                    thisExp.status = PAUSED
                    pauseExperiment(
                        thisExp=thisExp, 
                        win=win, 
                        timers=[globalClock], 
                    )
                    # once done pausing, restore running status
                    trials.status = STARTED
                thisExp.nextEntry()
                
            # completed nTrials repeats of 'trials'
            trials.status = FINISHED
            
            if thisSession is not None:
                # if running in a Session with a Liaison client, send data up to now
                thisSession.sendExperimentData()
            # mark thisCondition as finished
            if hasattr(thisCondition, 'status'):
                thisCondition.status = FINISHED
            # if awaiting a pause, pause now
            if conditions.status == PAUSED:
                thisExp.status = PAUSED
                pauseExperiment(
                    thisExp=thisExp, 
                    win=win, 
                    timers=[globalClock], 
                )
                # once done pausing, restore running status
                conditions.status = STARTED
            thisExp.nextEntry()
            
        # completed 1.0 repeats of 'conditions'
        conditions.status = FINISHED
        
        if thisSession is not None:
            # if running in a Session with a Liaison client, send data up to now
            thisSession.sendExperimentData()
        # mark thisBlock as finished
        if hasattr(thisBlock, 'status'):
            thisBlock.status = FINISHED
        # if awaiting a pause, pause now
        if blocks.status == PAUSED:
            thisExp.status = PAUSED
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[globalClock], 
            )
            # once done pausing, restore running status
            blocks.status = STARTED
    # completed 1.0 repeats of 'blocks'
    blocks.status = FINISHED
    
    
    # --- Prepare to start Routine "outro" ---
    # create an object to store info about Routine outro
    outro = data.Routine(
        name='outro',
        components=[outro_text],
    )
    outro.status = NOT_STARTED
    continueRoutine = True
    # update component parameters for each repeat
    outro_text.setText(task_texts['bye'])
    # store start times for outro
    outro.tStartRefresh = win.getFutureFlipTime(clock=globalClock)
    outro.tStart = globalClock.getTime(format='float')
    outro.status = STARTED
    outro.maxDuration = None
    # keep track of which components have finished
    outroComponents = outro.components
    for thisComponent in outro.components:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "outro" ---
    thisExp.currentRoutine = outro
    outro.forceEnded = routineForceEnded = not continueRoutine
    while continueRoutine and routineTimer.getTime() < 2.0:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *outro_text* updates
        
        # if outro_text is starting this frame...
        if outro_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            outro_text.frameNStart = frameN  # exact frame index
            outro_text.tStart = t  # local t and not account for scr refresh
            outro_text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(outro_text, 'tStartRefresh')  # time at next scr refresh
            # update status
            outro_text.status = STARTED
            outro_text.setAutoDraw(True)
        
        # if outro_text is active this frame...
        if outro_text.status == STARTED:
            # update params
            pass
        
        # if outro_text is stopping this frame...
        if outro_text.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > outro_text.tStartRefresh + 2-frameTolerance:
                # keep track of stop time/frame for later
                outro_text.tStop = t  # not accounting for scr refresh
                outro_text.tStopRefresh = tThisFlipGlobal  # on global time
                outro_text.frameNStop = frameN  # exact frame index
                # update status
                outro_text.status = FINISHED
                outro_text.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if defaultKeyboard.getKeys(keyList=["escape"]):
            thisExp.status = FINISHED
        if thisExp.status == FINISHED or endExpNow:
            endExperiment(thisExp, win=win)
            return
        # pause experiment here if requested
        if thisExp.status == PAUSED:
            pauseExperiment(
                thisExp=thisExp, 
                win=win, 
                timers=[routineTimer, globalClock], 
                currentRoutine=outro,
            )
            # skip the frame we paused on
            continue
        
        # has a Component requested the Routine to end?
        if not continueRoutine:
            outro.forceEnded = routineForceEnded = True
        # has the Routine been forcibly ended?
        if outro.forceEnded or routineForceEnded:
            break
        # has every Component finished?
        continueRoutine = False
        for thisComponent in outro.components:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "outro" ---
    for thisComponent in outro.components:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # store stop times for outro
    outro.tStop = globalClock.getTime(format='float')
    outro.tStopRefresh = tThisFlipGlobal
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if outro.maxDurationReached:
        routineTimer.addTime(-outro.maxDuration)
    elif outro.forceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-2.000000)
    thisExp.nextEntry()
    
    # mark experiment as finished
    endExperiment(thisExp, win=win)
    # end 'rush' mode
    core.rush(enable=False)


def saveData(thisExp):
    """
    Save data from this experiment
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    """
    filename = thisExp.dataFileName
    # these shouldn't be strictly necessary (should auto-save)
    thisExp.saveAsWideText(filename + '.csv', delim='semicolon')
    thisExp.saveAsPickle(filename)


def endExperiment(thisExp, win=None):
    """
    End this experiment, performing final shut down operations.
    
    This function does NOT close the window or end the Python process - use `quit` for this.
    
    Parameters
    ==========
    thisExp : psychopy.data.ExperimentHandler
        Handler object for this experiment, contains the data to save and information about 
        where to save it to.
    win : psychopy.visual.Window
        Window for this experiment.
    """
    if win is not None:
        # remove autodraw from all current components
        win.clearAutoDraw()
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed
        win.flip()
    # return console logger level to WARNING
    logging.console.setLevel(logging.WARNING)
    # mark experiment handler as finished
    thisExp.status = FINISHED
    # run any 'at exit' functions
    for fcn in runAtExit:
        fcn()
    logging.flush()


def quit(thisExp, win=None, thisSession=None):
    """
    Fully quit, closing the window and ending the Python process.
    
    Parameters
    ==========
    win : psychopy.visual.Window
        Window to close.
    thisSession : psychopy.session.Session or None
        Handle of the Session object this experiment is being run from, if any.
    """
    thisExp.abort()  # or data files will save again on exit
    # make sure everything is closed down
    if win is not None:
        # Flip one final time so any remaining win.callOnFlip() 
        # and win.timeOnFlip() tasks get executed before quitting
        win.flip()
        win.close()
    logging.flush()
    if thisSession is not None:
        thisSession.stop()
    # terminate Python process
    core.quit()


# if running this experiment as a script...
if __name__ == '__main__':
    # call all functions in order
    expInfo = showExpInfoDlg(expInfo=expInfo)
    thisExp = setupData(expInfo=expInfo)
    logFile = setupLogging(filename=thisExp.dataFileName)
    win = setupWindow(expInfo=expInfo)
    setupDevices(expInfo=expInfo, thisExp=thisExp, win=win)
    run(
        expInfo=expInfo, 
        thisExp=thisExp, 
        win=win,
        globalClock='float'
    )
    saveData(thisExp=thisExp)
    quit(thisExp=thisExp, win=win)
