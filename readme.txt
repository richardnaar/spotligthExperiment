Replication of the task in:

Müller, M. M., Malinowski, P., Gruber, T., & Hillyard, S. A. (2003). Sustained division of the attentional spotlight. Nature, 424(6946), 309-312. https://doi.org/10.1038/nature01812

# Experimental structure

The task is organized into blocks of 40 trials. Blocks are presented in a pseudorandom order, constrained so that the same experimental condition is not presented in more than two consecutive blocks.

At the beginning of each block, the starting response hand is selected randomly. The response hand is switched halfway through the block, after trial 20. An instruction screen is presented before each block and again before the halfway hand switch. The screen specifies:

* the locations to be attended (e.g., locations 1 and 3);
* the response hand to be used;
* whether a tone indicating the presence of a target pair will be presented;
* whether correct-response feedback will be provided; and
* the number of slowed-down trials at the beginning of the block.

# Training session

A separate training session can be run by setting trainingSession to 1 in the Info Dialog. The session lasts approximately one hour and should be completed on the same day as the main experiment, without EEG recording (EEG = 0).

The training session consists of 12 blocks. The first eight blocks provide progressively reduced auditory assistance:

During trials 1–20, a lower tone indicates the presence of a target pair, while a higher tone provides feedback following a correct response.
During trials 21–40, the target-pair tone is removed and only correct-response feedback is provided.

The first eight blocks also begin with five slowed-down trials. During these trials, each stimulus is presented for approximately 550 ms instead of approximately 183.3 ms. This extends the total trial duration from approximately 3.3 s to 9.9 s. Thus, the training session contains 40 slowed-down trials in total.

The final four blocks, comprising the last third of the training session, are intended to familiarize participants with the conditions of the main experiment. These blocks use the original stimulus timing and contain no target-pair tones, correct-response feedback, or slowed-down trials.

Make sure to:
- ... Configure the monitor settings in PsychoPy to match the required specifications, including screen size, viewing distance, and resolution. Click on the monitor icon to adjust these settings (by default, the experiment uses monitor settings saved under the name 'spotlight').
- ... To create new monitor settings, click on the monitor icon in PsychoPy and select New. If you already have a suitable monitor profile saved under a different name, you can update the experiment to use that by going to the Experiment Settings (ratchet icon) > Screen > Monitor.
- ... Set the monitor refresh rate to either 60 or 120 Hz and make the appropriate selection in the info dialog prior to the experiment.
- ... EEG:
	For the main EEG session, set trainingSession to 0 and EEG to 1.
	For the separate training session, set trainingSession to 1 and EEG to 0.
- ... Set EEG to 1 in the info dialog (it is off by default, but this can be changed by changing the order in the Experiment Info settings)
- ... Change the port address (the default can be modified in the Experiment info settings)
- ... Configure the audio output device in PsychoPy.Open the Device Manager (speaker icon in the top toolbar) and ensure that the device named Speakers corresponds to the actual physical speakers or headphones used in the experiment (e.g., Speakers (Realtek(R) Audio)). If necessary, add the correct device via Add device, remove unused devices, and keep the latency mode set to Shared low-latency (recommended).
- ... Set eyeTracker to 1 only when eye tracking is used. When enabled, verify that Experiment Settings (ratchet icon) > Data > Save hdf5 file is selected. Note: If the relevant eye-tracker plugin is not installed, the Eyetracker device menu may contain only None and MouseGaze. The specific eye tracker option should become available after the plugin has been installed and PsychoPy has been restarted.

# Eye-tracking setup

Support for physical eye trackers is provided through device-specific PsychoPy plugins and may not be included in the default installation. To configure a Tobii eye tracker (similar to other trackers such as GazePoint):

Ensure that the required Tobii software and drivers are installed and verify that the tracker is detected independently of PsychoPy.
In PsychoPy, open Tools > Plugin/packages manager...
Select the Plugins tab, search for Tobii Eyetracker Support, and install psychopy-eyetracker-tobii.
Restart PsychoPy after installation.
Open Experiment Settings (ratchet icon) > Eyetracking and select Tobii Technology from the Eyetracker device menu. Configure the appropriate tracker model and any additional device-specific settings.
Open Experiment Settings > Data and enable Save hdf5 file to store the eye-tracking data.
Set eyeTracker to 1 in the Info Dialog when eye tracking is used. Set it to 0 for runs without eye tracking.

Installing psychopy-eyetracker-tobii should also install the required tobii_research dependency automatically. The plugin is then loaded when PsychoPy restarts. The HDF5 option is intended for complex data such as eye-tracking recordings.

---

To ensure that on-screen objects correspond to the intended degrees of visual angle:

- ... Set the viewing distance to 57 cm for debuging. At this distance, degrees and centimeters should match, allowing you to measure objects directly on the screen (e.g., the box should measure 2.5 x 3.2 cm).

---

# Adding and Managing Translations

This experiment supports multiple languages through external translation files located in the translations/ directory. Each translation file follows the naming convention:

For example:

translation_ENG.xlsx — English
translation_EST.xlsx — Estonian

To add a new language:

Create a new .xlsx file in the translations/ directory, following the existing format.

Name the file as translation_<LANGUAGE_CODE>.xlsx, replacing <LANGUAGE_CODE> with a short, clear identifier for the language (e.g., FR for French).

Ensure the file structure matches the existing translation files:

A column named key for reference terms.
A column named translated_message for the corresponding translations.
Select the new language from the Info Dialog at the start of the experiment.

Note: The experiment will use the file corresponding to the language selected in the translation field in the Info Dialog.

---

TestRun Mode:
- When the "testRun" option is enabled in the Info Dialog, the experiment simulates keyboard responses. For example, simulated responses (dummy responses) are triggered at predetermined reaction times from a list (RTs), so that you can test the flow of the experiment without requiring manual key presses.
- Note that instructions within the Test Mode are presented for only a brief duration and continues automatically after the allotted time.

---

# DATA FILES AND VARIABLES

The experiment produces several output files per run. The primary file used for analysis is the .csv file, although the same information (and additional metadata) is also stored in the .psydat file. The .log file is a plain-text file that provides a chronological record of events during the run (e.g., stimulus presentation, timing information, and warnings).

The .csv file contains, among others, the following key columns:

- attendCond: experimental condition
- isTraining: indicates whether a trial belongs to the training phase (1 = training)
- targetPair: indicates whether a target was present (TRUE / FALSE)
- accuracy: response accuracy
- RT: response classification and reaction time

For non-training performance, the target detection rate (TDR) can be computed from the subset of trials where isTraining == 0 and targetPair == TRUE, using the accuracy column.

The RT column encodes response outcomes as follows: hits are represented by a decimal value corresponding to the reaction time, misses are labeled "miss" and count as incorrect responses, and false alarms are labeled "false alarm". While false alarms are counted as incorrect responses, they are not included in the TDR calculation when considering only trials in which a target was present (i.e., targetPair == TRUE).

---

This experiment was created using PsychoPy3 Experiment Builder (v2026.1.1)

Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
    PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
    https://doi.org/10.3758/s13428-018-01193-y