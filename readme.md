# Müller et al. (2003) split-attention task replication

This experiment replicates the task described in:

> Müller, M. M., Malinowski, P., Gruber, T., & Hillyard, S. A. (2003). Sustained division of the attentional spotlight. *Nature, 424*(6946), 309–312. [https://doi.org/10.1038/nature01812](https://doi.org/10.1038/nature01812)

## Experimental structure

The task is organized into blocks of 40 trials. Blocks are presented in a pseudorandom order, constrained so that the same experimental condition is not presented in more than two consecutive blocks.

At the beginning of each block, the starting response hand is selected at random. The response hand is switched halfway through the block, after trial 20. An instruction screen is presented before each block and again immediately before the halfway hand switch. The screen specifies:

- the locations to be attended (e.g., locations 1 and 3);
- the response hand to be used;
- whether a tone indicating the presence of a target pair will be presented;
- whether correct-response feedback will be provided; and
- the number of slowed-down trials at the beginning of the block.

## Training session

A separate training session can be run by setting `trainingSession` to `1` in the Experiment Info dialog. The session lasts approximately one hour and should be completed on the same day as the main experiment, without EEG recording (`EEG = 0`).

The training session consists of 12 blocks. The first eight blocks include auditory assistance, which is reduced halfway through each block:

- During trials 1–20, a lower tone indicates the presence of a target pair, while a higher tone provides feedback following a correct response.
- During trials 21–40, the target-pair tone is removed and only correct-response feedback is provided.

The first eight blocks also begin with five slowed-down trials. During these trials, each stimulus is presented for approximately 550 ms instead of approximately 183.3 ms. This extends the total trial duration from approximately 3.3 s to 9.9 s. Thus, the training session contains 40 slowed-down trials in total.

The final four blocks, comprising the last third of the training session, are intended to familiarize participants with the conditions of the main experiment. These blocks use the original stimulus timing and contain no target-pair tones, correct-response feedback, or slowed-down trials.

## Experiment setup

### Monitor configuration

- Configure the monitor settings in PsychoPy to match the required screen size, viewing distance, and resolution. Click the monitor icon to adjust these settings. By default, the experiment uses the monitor profile named `spotlight`.
- To create a monitor profile, click the monitor icon and select **New**. To use an existing profile saved under another name, open **Experiment Settings (ratchet icon) > Screen > Monitor** and select it.
- Set the monitor refresh rate to either 60 or 120 Hz and select the corresponding `refreshRate` value in the Experiment Info dialog before starting the experiment.

### Session and EEG configuration

- For the main EEG session, set `trainingSession` to `0` and `EEG` to `1`.
- For the separate training session, set `trainingSession` to `1` and `EEG` to `0`.
- Verify the EEG port address before recording. The default can be changed in the Experiment Info settings.

> **Note:** `EEG` is set to `0` by default. This default can be changed by reordering the choices in the Experiment Info settings.

### Audio configuration

Open the PsychoPy Device Manager using the speaker icon in the top toolbar. Ensure that the device named `Speakers` corresponds to the physical speakers or headphones used in the experiment (e.g., `Speakers (Realtek(R) Audio)`). If necessary, select **Add device**, add the correct output device, and remove unused devices. Keep the latency mode set to **Shared low-latency (recommended)**.

### Eye-tracking setup

Support for physical eye trackers is provided through device-specific PsychoPy plugins and may not be included in the default installation. The following steps apply to a Tobii eye tracker. For another supported tracker, such as Gazepoint, install and select the corresponding plugin instead.

1. Ensure that the required Tobii software and drivers are installed, and verify that the tracker is detected independently of PsychoPy.
2. In PsychoPy, open **Tools > Plugin/packages manager...**
3. Select the **Plugins** tab, search for **Tobii Eyetracker Support**, and install `psychopy-eyetracker-tobii`.
4. **Restart PsychoPy after installing the plugin.**
5. Open **Experiment Settings (ratchet icon) > Eyetracking** and select **Tobii Technology** from the **Eyetracker device** menu. Configure the appropriate tracker model and any additional device-specific settings.
6. Open **Experiment Settings > Data** and enable **Save hdf5 file** to store the eye-tracking data.
7. Set `eyeTracker` to `1` in the Experiment Info dialog when eye tracking is used. Set it to `0` for runs without eye tracking.

If the relevant plugin is not installed, the **Eyetracker device** menu may contain only **None** and **MouseGaze**. The physical eye-tracker option should become available after the plugin has been installed and PsychoPy has been restarted.

Installing `psychopy-eyetracker-tobii` should also install its required `tobii_research` dependency automatically. The plugin is loaded when PsychoPy restarts. The HDF5 option is intended for complex data such as eye-tracking recordings.

For additional information, see PsychoPy's [eye-tracking documentation](https://psychopy.org/hardware/eyeTracking.html) and the [`psychopy-eyetracker-tobii` plugin documentation](https://github.com/psychopy/psychopy-eyetracker-tobii).

### Visual-angle verification

For debugging, set the viewing distance to 57 cm. At this distance, 1° of visual angle corresponds to approximately 1 cm for small objects, allowing the dimensions of on-screen objects to be measured directly. For example, a 2.5° × 3.2° box should measure approximately 2.5 × 3.2 cm on the screen.

## Adding and managing translations

This experiment supports multiple languages through external translation files in the `translations/` directory. Translation files follow the naming convention `translation_<LANGUAGE_CODE>.xlsx`, for example:

- `translation_ENG.xlsx` — English
- `translation_EST.xlsx` — Estonian

To add a language:

1. Create an `.xlsx` file in the `translations/` directory, following the format of the existing translation files.
2. Name it `translation_<LANGUAGE_CODE>.xlsx`, replacing `<LANGUAGE_CODE>` with a short, clear identifier (e.g., `translation_FR.xlsx` for French).
3. Include the following columns:
   - `key` for the reference terms;
   - `translated_message` for the corresponding translations.
4. Add the language code to the choices for `translation` in the Experiment Info settings.
5. Select the language from the Experiment Info dialog when starting the experiment.

> **Note:** The value selected in the `translation` field must correspond to the language code in the translation filename.

## Test-run mode

- When `testRun` is set to `1` in the Experiment Info dialog, the experiment simulates keyboard responses. Dummy responses are triggered at predetermined reaction times from the `RTs` list, allowing the experimental flow to be tested without manual keypresses.
- Instructions in test-run mode are presented only briefly and continue automatically after the allotted time.

## Data files and variables

The experiment produces several output files per run. The primary file used for analysis is the `.csv` file, although the same information and additional metadata are also stored in the `.psydat` file. The `.log` file is a plain-text chronological record of events during the run, including stimulus presentation, timing information, and warnings.

The `.csv` file contains, among others, the following key columns:

| Column | Description |
| --- | --- |
| `attendCond` | Experimental condition |
| `isTraining` | Whether the trial belongs to the training phase (`1` = training) |
| `targetPair` | Whether a target was present (`TRUE` or `FALSE`) |
| `accuracy` | Response accuracy |
| `RT` | Response classification and reaction time |

For non-training performance, the target detection rate (TDR) can be computed from trials for which `isTraining == 0` and `targetPair == TRUE`, using the `accuracy` column.

The `RT` column encodes response outcomes as follows:

| `RT` value | Meaning |
| --- | --- |
| Decimal value | Hit; the value is the reaction time |
| `"miss"` | Miss; counted as an incorrect response |
| `"false alarm"` | False alarm; counted as an incorrect response |

False alarms are not included in the TDR because the TDR is calculated only from trials in which a target was present (`targetPair == TRUE`).

## Software and citation

This experiment was created using PsychoPy Builder v2026.1.1.

> Peirce, J. W., Gray, J. R., Simpson, S., MacAskill, M., Höchenberger, R., Sogo, H., Kastman, E., & Lindeløv, J. K. (2019). PsychoPy2: Experiments in behavior made easy. *Behavior Research Methods, 51*, 195–203. [https://doi.org/10.3758/s13428-018-01193-y](https://doi.org/10.3758/s13428-018-01193-y)