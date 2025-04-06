clear all

mainPath    =   'C:\Users\pcadmin\Documents\dok\TARU\PROJEKTID\EEGManyLabs\Task\spotlight_replication\analysis\';
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

dati = 2;

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

elec2plot = 73; %erg1 == light sensor

condi = 2;
condIndx = find(contains({EEG.event.resp}, config.resp{condi}));
onset = [EEG.event(condIndx).latency];

% juku = [];
% for indx = 1:length(onset)-1
%     juku = [juku; EEG.event(condIndx(indx+1)).latency-EEG.event(condIndx(indx)).latency];
% end
% 
% juku1 = juku/EEG.srate;


meanComplexFFT = [];
for eventi = 1:length(onset)
    
    transient = 0; % in samples
    eventDur = 3.1*EEG.srate;
    current_data = EEG.data(elec2plot, round(onset(eventi)):round((onset(eventi)+(eventDur))) );
    fft_dat = fft(current_data);
    % meanComplexFFT(:,eventi) = fft_dat/length(current_data); % coherent
    meanComplexFFT(:,eventi) = abs(fft_dat/length(current_data));
end

% juku = abs(mean(meanComplexFFT,2)); % coherent
juku = mean(meanComplexFFT,2);

hz = 0:1/3.1:EEG.srate/2;
stem(hz(2:500), juku(2:500))

eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);


