clear all

mainPath    =   'D:\Müller2003Data\';
% mainPath    =   'C:\Users\pcadmin\Documents\dok\TARU\PROJEKTID\EEGManyLabs\Task\spotlight_replication\analysis\';

locfile     =   'C:\Users\pcadmin\Documents\dok\TARU\eeglab2023.0\64_4EOG.ced';                                     % dir of the electrode location file 

pSensorPath = [mainPath, 'photosensor\'];
eegPath = [mainPath, 'EEG\'];

addpath(pSensorPath)
addpath(mainPath)


eeglab(), close

%% trig

% Load the configuration
run('trigger_config_myller_Tartu.m');

% Check if the config variable exists
if ~exist('config', 'var')
    error('Configuration file did not load correctly or "config" variable is missing.');
end

%% Light sensor

impdir = pSensorPath;
implist = dir([impdir, '*.bdf']);

dati = 1; %  

[~, subName, ~] = fileparts(implist(dati).name);
    
ALLEEG = []; EEG = []; CURRENTSET = [];                                                                        % initialize EEGLAB data structures for the new round
        
EEG = pop_biosig([impdir, implist(dati).name]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, ...
    'setname', implist(dati).name); 
    
    
EEG = gtrigf(EEG,config, 0, 'stim');

triggerCounts = countTriggers(EEG, config);


% Down-sample
if EEG.srate ~= 500
    EEG = pop_resample(EEG, 500);
end

elec2plot = 73; % 73 erg1 == light sensor

condi = 2;
condIndx = find(contains({EEG.event.resp}, config.resp{condi}));
onset = [EEG.event(condIndx).latency];


meanComplexFFT = [];
for eventi = 1:length(onset)
    transient = 0; % in samples
    eventDur = 3.1*EEG.srate;
    current_data = EEG.data(elec2plot, round(onset(eventi)):round((onset(eventi)+(eventDur-1))) );
    fft_dat = fft(current_data);
    % meanComplexFFT(:,eventi) = fft_dat/length(current_data); % coherent
    meanComplexFFT(:,eventi) = abs(fft_dat)/length(current_data)*2;
end

% juku = abs(mean(meanComplexFFT,2)); % coherent
juku = mean(meanComplexFFT,2);

hz = (0:length(juku)-1) * (EEG.srate / length(juku));
% hz = 0:1/3.1:EEG.srate/2;
stem(hz(2:500), juku(2:500))
plot(hz(2:100), juku(2:100))

eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

%% This is a raw version, need to be cleaned and verified

impdir = eegPath; implist = dir([impdir, '*.bdf']);

subi = 1

subName = implist(subi).name;

fprintf('loading participant: %s \n', subName);                                                                     % print the file name to the Command Window

ALLEEG = []; EEG = []; CURRENTSET = [];    

EEG = pop_biosig([impdir, implist(subi).name]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, ...
    'setname', implist(subi).name);   


%% Triggers

EEG = gtrigf(EEG,config, 0, 'stim');

triggerCounts = countTriggers(EEG, config);
%% Frequencies based on the refresh rate

monitor_refresh = 119.89; % 60.8
multip = round(monitor_refresh/60); 
nFrames = [7,5,4,3]*multip;

stim_frex = 1000./(1000/monitor_refresh*nFrames);
symbol_frex = 1000./(1000/monitor_refresh*11*multip);


%%

exg1 = find(strcmp('EXG1', {EEG.chanlocs.labels}));

% compute mastoid reference
 
% EEG = pop_reref(EEG, [exg1+4 exg1+5]);
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% remove unnecessary
EEG.data = EEG.data([1:64 exg1:exg1+3],:); EEG.nbchan = 64; EEG.subject = implist(subi).name;          % remove unnecessary channels (64 head channels, 4 eyes, additional EMG channel)
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

EEG = pop_chanedit(EEG, 'load',{locfile 'filetype' 'autodetect'});                                     % edit channel locations structure of an EEGLAB dataset, EEG.chanlocs. 
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 


%% EOG electrodes
eogs = {'VEOG', 'HEOG'}; 

eogi = find(contains({EEG.chanlocs.labels}, eogs(1:2)));
EOG = EEG.data(eogi,:); EOGlocs = EEG.chanlocs(eogi);

EEG.data = EEG.data(1:end-length(eogi), :); EEG.chanlocs = EEG.chanlocs(1:end-length(eogi)); EEG.nbchan = 64;

%% Not implemented in the original
% Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal. 
EEG = clean_flatlines(EEG, 5);

% see EEG.chaninfo.removedchans 
[EEG, rej] = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',8,'norm','on','measure','kurt'); 
    
%% compute average reference
% eegplot(EEG.data, 'events', EEG.event) 
EEG_ref = mean(EEG.data(1:end,:),1); % 1:64

% not adding EOG channel back at the moment
% EEG.data = [EEG.data; EOG]; EEG.chanlocs(end+1) = EOGlocs; EEG.nbchan = EEG.nbchan + 1;

EEG.data = bsxfun(@minus, EEG.data(1:end,:), squeeze(EEG_ref)); % average reference (excluding EOG from the mean)

fprintf(['Averge reference (excluding EOG from the mean) \n'])


%% DOWN-SAMPLE 
if EEG.srate == 500
else
    EEG = pop_resample(EEG, 500);
end

%% filter 
lowcutoff = 1;
revfilt = 0; % [0|1] reverse filter (i.e. bandpass filter to notch filter). {default 0}
plotfreqz = 0;
minphase = false; % defalut

[EEG, com, b] = pop_eegfiltnew(EEG, lowcutoff, 125, [], revfilt, [], plotfreqz, minphase);

%%
% electrode = 'PO8';
pool = {'PO3','PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'O1', 'O2', 'I1', 'I2', 'Oz', 'Iz', 'POz'};
% pool = {'PO3','PO4','PO7','PO8','O1','O2'};

elec2plot = find(ismember({EEG.chanlocs.labels}, pool));

condi = 2;
condIndx = find(contains({EEG.event.resp}, config.resp{condi}));
onset = [EEG.event(condIndx).latency];

startInS = 0.526;
dur = 3.1;
meanComplexFFT = [];
for eventi = 1:length(onset)
    
    transient = startInS*EEG.srate; % in samples
    eventDur = dur*EEG.srate;
    % current_data = EEG.data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) );
    current_data = mean(EEG.data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) ),1); % avg over electrodes
    fft_dat = fft(current_data);
    meanComplexFFT(:,eventi) = fft_dat/length(current_data); % coherent
    % meanComplexFFT(:,eventi) = abs(fft_dat)/length(current_data)*2;
end


juku = abs(mean(meanComplexFFT,2)); % lose coherence
% juku = mean(meanComplexFFT,2);

hz = 0:1/(dur-startInS):EEG.srate/2;
% hz = (0:length(juku)-1) * (EEG.srate / length(juku));

stem(hz(2:500), juku(2:500))
plot(hz(2:100), juku(2:100))

% eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

skipbins = 1; 
numbins = 1 + skipbins;

% Computing SNR
for hzi = numbins + 1:length(juku) - numbins - 1
    numer = juku(hzi);
    denom = rms(juku([hzi-numbins:hzi-skipbins, hzi+skipbins:hzi+numbins]));
    snrE(hzi) = numer./denom;
end 

plot(hz(2:100), snrE(2:100))

