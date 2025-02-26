Replication of the task in:

Müller, M. M., Malinowski, P., Gruber, T., & Hillyard, S. A. (2003). Sustained division of the attentional spotlight. Nature, 424(6946), 309-312. https://doi.org/10.1038/nature01812


Make sure to:
- ... Set the monitor refresh rate to either 60 or 120 Hz and make the appropriate selection in the info dialog prior to the experiment.
- ... Configure the monitor settings in PsychoPy to match the required specifications, including screen size, viewing distance, and resolution. Click on the monitor icon to adjust these settings (by default, the experiment uses monitor settings saved under the name 'spotlight').
- ... Set EEG to 1 in the info dialog (it is off by default, but this can be changed by changing the order in the Experiment Info settings)
- ... Change the port address (the default can be modified in the Experiment info settings)

---

To ensure that on-screen objects correspond to the intended degrees of visual angle:

- ... Set the viewing distance to 57 cm for debuging. At this distance, degrees and centimeters should match, allowing you to measure objects directly on the screen (e.g., the box should measure 2.5 x 3.2 cm)."

---

Adding and Managing Translations

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

This experiment was created using PsychoPy3 Experiment Builder (v2024.2.1)

Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
    PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
    https://doi.org/10.3758/s13428-018-01193-y