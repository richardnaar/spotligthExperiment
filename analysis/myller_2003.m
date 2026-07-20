clear all

mainPath    =   'D:\Müller2003Data\';
codePath    =   'C:\Users\pcadmin\Documents\dok\TARU\PROJEKTID\EEGManyLabs\Task\spotlight_replication\analysis\';

locfile     =   'C:\Users\pcadmin\Documents\dok\TARU\eeglab2023.0\64_4EOG.ced';                                     % dir of the electrode location file 

pSensorPath = [mainPath, 'photosensor\'];
eegPath = [mainPath, 'EEG\'];

addpath(pSensorPath)
addpath(mainPath)

cd(codePath)

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
    
    
EEG = gtrigf(EEG,config, 0, 'start');

triggerCounts = countTriggers(EEG, config);


% Down-sample
if EEG.srate ~= 500
    EEG = pop_resample(EEG, 500);
end

elec2plot = 73; % 73 erg1 == light sensor

condi = 1; starti = 1;
cond = contains({EEG.event.cond}, config.cond{condi});
start = contains({EEG.event.position}, config.position{starti});
condIndx = find(cond.*start);
onset = [EEG.event(condIndx).latency];


meanComplexFFT = [];
for eventi = 1:length(onset)
    transient = EEG.srate*0.526; % in samples
    eventDur = 3.1*EEG.srate;
    current_data = EEG.data(elec2plot, round(onset(eventi)):round((onset(eventi)+(eventDur-1))) );
    fft_dat = fft(current_data);
    meanComplexFFT(:,eventi) = fft_dat/length(current_data); % coherent
    % meanComplexFFT(:,eventi) = abs(fft_dat)/length(current_data)*2;
end

juku = abs(mean(meanComplexFFT,2)); % coherent
% juku = mean(meanComplexFFT,2);

hz = (0:length(juku)-1) * (EEG.srate / length(juku));
% hz = 0:1/3.1:EEG.srate/2;
% stem(hz(2:100), juku(2:100))
plot(hz(2:100), juku(2:100))

% eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

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

EEG = gtrigf(EEG,config, 0, 'start');

triggerCounts = countTriggers(EEG, config);
%% Frequencies based on the refresh rate

monitor_refresh = 119.89; % 60.8
multip = round(monitor_refresh/60); 
% multip = monitor_refresh/60; 
nFrames = [7,5,4,3]*multip;

stim_frex = 1000./(1000/monitor_refresh*nFrames);
symbol_frex = 1000./(1000/monitor_refresh*11*multip);


%%

exg1 = find(strcmp('EXG1', {EEG.chanlocs.labels}));

% compute mastoid reference. here?
 
EEG = pop_reref(EEG, [exg1+4 exg1+5]);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

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
    
%% compute average reference *not in this project
% % eegplot(EEG.data, 'events', EEG.event) 
% EEG_ref = mean(EEG.data(1:end,:),1); % 1:64
% 
% % not adding EOG channel back at the moment
% % EEG.data = [EEG.data; EOG]; EEG.chanlocs(end+1) = EOGlocs; EEG.nbchan = EEG.nbchan + 1;
% 
% EEG.data = bsxfun(@minus, EEG.data(1:end,:), squeeze(EEG_ref)); % average reference (excluding EOG from the mean)
% 
% fprintf(['Averge reference (excluding EOG from the mean) \n'])


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
% pool = {'PO3','PO4', 'PO7', 'PO8', 'O1', 'O2'};
% pool = {'PO3','PO4','PO7','PO8','O1','O2'};
pool =  {'PO8', 'PO7'};

elec2plot = find(ismember({EEG.chanlocs.labels}, pool));

% condi = 4;
for condi = 1:4
start_idx = contains({EEG.event.position}, config.position{1});
cond_idx = contains({EEG.event.cond}, config.cond{condi});
condIndx = find(start_idx.*cond_idx);

% condIndx = find(contains({EEG.event.position}, config.position{1}));
onset = [EEG.event(condIndx).latency];

startInS = 0.526;
dur = 3.1;
meanComplexFFT = [];
for eventi = 1:length(onset)
    transient = startInS*EEG.srate; % in samples
    eventDur = dur*EEG.srate;
    current_data = mean(EEG.data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) ),1); % avg over electrodes
    fft_dat = fft(current_data);
    meanComplexFFT(:,eventi) = fft_dat; % coherent
end


juku = abs(mean(meanComplexFFT,2))/size(meanComplexFFT,1); % lose coherence
% juku = mean(meanComplexFFT,2);

hz = 0:1/(dur-startInS):EEG.srate/2;
% hz = (0:length(juku)-1) * (EEG.srate / length(juku));

% stem(hz(2:500), juku(2:500))
plot(hz(2:100), juku(2:100))

% eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

skipbins = 1; 
numbins = 3 + skipbins;

% Computing SNR
for hzi = numbins + 1:length(juku) - numbins - 1
    numer = juku(hzi);
    denom = rms(juku([hzi-numbins:hzi-skipbins, hzi+skipbins:hzi+numbins]));
    snrE(hzi) = numer./denom;
end 

plot(hz(2:100), snrE(2:100))

%% complex demod 
params.targetFreqs    = stim_frex;     % Hz; SSVEP frequencies of interest

params.epochPre       = -0.5;                    % seconds relative to event time (e.g. sync or flicker onset)
params.epochPost      =  3.5;                    % seconds after event

% Taper (square–cosine ramp) as fraction of window length at each edge
% e.g. 0.1 = 10% of the samples at each edge tapered
params.taperFraction  = 0.1;      

% Low-pass for envelope (how “slowly” amplitude can vary)
params.lpCutoffHz     = 1.0;                     % Hz (0.5–2 Hz is typical for SSVEP envelope)
params.lpOrder        = 3;                       % Butterworth order (3–6 is usually fine)

params.analysisStart  = 0.526;                   % s after flicker onset
params.analysisEnd    = 3.1;                     % 3.026;  % s after flicker onset

fs        = EEG.srate;
preSamp   = round(params.epochPre  * fs);  % negative
postSamp  = round(params.epochPost * fs);  % positive
sampleVec = preSamp:postSamp;
nSamp     = numel(sampleVec);

% time axis relative to event (e.g. flicker onset or sync)
tEpoch = sampleVec / fs;   % e.g. -0.5 … +3.5
nTrials = numel(onset);

% Allocate: ROI-averaged epochs (time x trial)
roiEpochs = nan(nSamp, nTrials);

% ------- extract epochs and average across ROI electrodes -------
for ti = 1:nTrials
    zeroIdx = round(onset(ti));  % event sample
    segIdx  = zeroIdx + sampleVec;

    % guard against running off the data edges
    if segIdx(1) < 1 || segIdx(end) > size(EEG.data,2)
        continue;
    end

    % mean over ROI electrodes, store as time x trial
    roiEpochs(:, ti) = squeeze(mean(EEG.data(elec2plot, segIdx), 1));
end

% average waveform across trials (this is what the paper describes)
% result: 1 x nSamp (or nSamp x 1)
avgWaveform = nanmean(roiEpochs, 2)';  % make it row vector

% ------- complex demodulation of the *averaged waveform* -------
targetFreqs = params.targetFreqs;
nFreqs      = numel(targetFreqs);

demAmpFreq   = nan(nFreqs, nSamp);
demPhaseFreq = nan(nFreqs, nSamp);

for fi = 1:nFreqs
    [demAmp, demPhase] = complex_demod_segment( ...
        avgWaveform, fs, targetFreqs(fi), ...
        params.lpCutoffHz, params.lpOrder, params.taperFraction);

    demAmpFreq(fi, :)   = demAmp;
    demPhaseFreq(fi, :) = demPhase;
end

% ------- pick the analysis window (exclude first 500 ms, etc.) -------

% Here, tEpoch is relative to the event. If your event is at flicker onset,
% then tEpoch already matches the paper’s timing.
% If your event is at "sync" = 26 ms after flicker, you can either:
% (a) adjust params.analysisStart/End as noted above, or
% (b) build a time axis explicitly relative to flicker and use that.


plot(demAmpFreq')
pause()
%%

tFlicker = tEpoch;  % adjust this if event != flicker; see comment in params

mask = tFlicker >= params.analysisStart & tFlicker <= params.analysisEnd;

% mean SSVEP amplitude over analysis interval, per frequency
meanAmp_scalar = mean(demAmpFreq(:, mask), 2);   % nFreqs x 1

% store results
grandDemod(subi).cond(condi).tEpoch           = tEpoch;
grandDemod(subi).cond(condi).freqs           = targetFreqs;
grandDemod(subi).cond(condi).avgWaveform     = avgWaveform;
grandDemod(subi).cond(condi).demAmpFreq      = demAmpFreq;
grandDemod(subi).cond(condi).demPhaseFreq    = demPhaseFreq;
grandDemod(subi).cond(condi).meanAmp_scalar  = meanAmp_scalar;
end
%% fi 2 and 3 at conditions 
% attended adjacent 1+2 for 3 and 3+4 for 2; 
% attended separate 1+3 for 3 and 2+4 for 2; 
% unattended separate 1+3 for 2 and 2+4 for 3;
% unattended oppposite 1+2 for 3 and 3+4 for 2;

% save([output, 'grandDemod-', date, '.mat'], 'grandDemod', '-v7.3');
% outlist = dir([output, '*.mat']);
% load([output, outlist(1).name])
%%

% ----- collect amplitudes from grandDemod into an array -----
nSub  = numel(grandDemod);
nCond = numel(grandDemod(1).cond);
nFreq = numel(grandDemod(1).cond(1).meanAmp_scalar);

allAmp = nan(nSub, nCond, nFreq);   % subj × cond × freq

for subi = 1:nSub
    for condi = 1:nCond
        allAmp(subi,condi,:) = grandDemod(subi).cond(condi).meanAmp_scalar(:);
    end
end

% ----- mean and SEM across subjects -----
avg_dat = squeeze(mean(allAmp,1));          % cond × freq
avg_dat = avg_dat.';                        % freq × cond (to match your old code)

N   = nSub;
err = squeeze(std(allAmp,[],1)).'/sqrt(N);  % freq × cond  (SEM)

freqs = grandDemod(1).cond(1).freqs; 




%% plot

% Define colors: first 3 blue, last 3 red
% colors = [repmat([0 0 1], 2, 1); repmat([1 0 0], 2, 1)];

colors = [repmat([0 0 0], 4, 1)];

defCondOrder = {'1+2', '1+3', '2+4', '3+4'};
conOrder2 = [[1,2,3,4]; [1,3,2,4]; [4,2,3,1]; [4,3,2,1]];   % default condition order: '1+2', '1+3', '2+4', '3+4', frex order from lowest to higest (e.g. 8.7, 12.2, 15.2, 20.3)
forder = [3,1,4,2];
figure;
for f = 1:size(avg_dat, 1)    % loop over frequencies
    
    f2show = forder(f);
    data = avg_dat(f2show, conOrder2(f,:));
    e    = err(f2show,  conOrder2(f,:));

    b = bar(data, 'FaceColor', 'flat');
    hold on
    for c = 1:length(conOrder2(f,:))
        b.CData(c,:) = colors(c,:);
    end
    
    ylabel('Amplitude (\mu V)'); % 
    % xlabel('Condition');
    xticks(1:6);
    xticklabels(conOrder2(f,:));
    % ylim([0 max(avg_dat(:))*1.1]);
    
    % error bars (using your errorb function)
    % errorb(1:4, data, e, 'top');

    % optional: add frequency in the title
    % freqs(f) assumed from above
    title(sprintf('Position %.0f (%.1f Hz)', [f, freqs(f2show)]));
    xticks(1:4)
    xticklabels(defCondOrder(conOrder2(f,:)))
    % defCondOrder{conOrder2(f,:)}

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','fontname'),'fontname', 'Arial')

    pause();   % step through frequencies
    hold off
end

%%

% allAmp(subi, condi, frexi)

defCondOrder = {'1+2', '1+3', '2+4', '3+4'};

condIdx = @(name) find(strcmp(defCondOrder, name));

c12 = condIdx('1+2');
c13 = condIdx('1+3');
c24 = condIdx('2+4');
c34 = condIdx('3+4');

% Target frequencies / positions
pos2_frex = 1;   % e.g. 8.7 Hz, position 2
pos3_frex = 4;   % e.g. 20.3 Hz, position 3

% Max across conditions, separately for each subject and frequency
% Result: nSub x 1 x nFrex
subFreqMax = max(allAmp, [], 2, 'omitnan');

subFreqMax(subFreqMax == 0) = NaN;

% Normalized amplitudes
allAmpNorm = allAmp ./ subFreqMax;

%%

% -----------------------------
% Attended adjacent
% position 2 attended in 1+2
% position 3 attended in 3+4
% -----------------------------
attendedAdjacent = mean([ ...
    allAmpNorm(:, c12, pos2_frex), ...
    allAmpNorm(:, c34, pos3_frex)], ...
    2, 'omitnan');

% -----------------------------
% Attended separated
% position 2 attended in 2+4
% position 3 attended in 1+3
% -----------------------------
attendedSeparated = mean([ ...
    allAmpNorm(:, c24, pos2_frex), ...
    allAmpNorm(:, c13, pos3_frex)], ...
    2, 'omitnan');

% -----------------------------
% Unattended intermediate
% position 2 ignored in 1+3
% position 3 ignored in 2+4
% -----------------------------
unattendedIntermediate = mean([ ...
    allAmpNorm(:, c13, pos2_frex), ...
    allAmpNorm(:, c24, pos3_frex)], ...
    2, 'omitnan');

% -----------------------------
% Unattended opposite
% position 2 ignored in 3+4
% position 3 ignored in 1+2
% -----------------------------
unattendedOpposite = mean([ ...
    allAmpNorm(:, c34, pos2_frex), ...
    allAmpNorm(:, c12, pos3_frex)], ...
    2, 'omitnan');

plotVals = [ ...
    attendedAdjacent, ...
    attendedSeparated, ...
    unattendedIntermediate, ...
    unattendedOpposite];

%%

catLabels = { ...
    'Attended adjacent', ...
    'Attended separated', ...
    'Unattended intermediate', ...
    'Unattended opposite'};

% Put each word on a separate line
catLabels = strrep(catLabels, ' ', sprintf('\n'));

groupMean = mean(plotVals, 1, 'omitnan');
groupSEM  = std(plotVals, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(plotVals), 1));

figure;
bar(groupMean, 'k');
hold on;

% errorbar(1:4, groupMean, groupSEM, ...
%     'k', ...
%     'LineStyle', 'none', ...
%     'LineWidth', 1.5);

ax = gca;

% Keep tick positions, but remove MATLAB's default tick labels
ax.XTick = 1:4;
ax.XTickLabel = [];

ylabel('Normalized amplitude  (\mu V)');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(findall(gcf,'-property','fontname'),'fontname', 'Arial')
ylim([0 1.05]);
box off;

% Make extra space below the axis for custom labels
ax.Position(2) = ax.Position(2) + 0.10;
ax.Position(4) = ax.Position(4) - 0.10;

% Add custom multiline x-axis labels
yl = ylim(ax);
yText = yl(1) - 0.02 * range(yl);

for i = 1:numel(catLabels)
    text(i, yText, catLabels{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', ax.FontSize, ...
        'Clipping', 'off');
end

%%

% Attended adjacent vs attended separated
[~, p_adj_sep, ~, stats_adj_sep] = ttest(attendedAdjacent, attendedSeparated);

% Attended separated vs unattended intermediate
[~, p_sep_int, ~, stats_sep_int] = ttest(attendedSeparated, unattendedIntermediate);

% Unattended intermediate vs unattended opposite
[~, p_int_opp, ~, stats_int_opp] = ttest(unattendedIntermediate, unattendedOpposite);

fprintf('Attended adjacent vs attended separated: t(%d)=%.3f, p=%.5f\n', ...
    stats_adj_sep.df, stats_adj_sep.tstat, p_adj_sep);

fprintf('Attended separated vs unattended intermediate: t(%d)=%.3f, p=%.5f\n', ...
    stats_sep_int.df, stats_sep_int.tstat, p_sep_int);

fprintf('Unattended intermediate vs unattended opposite: t(%d)=%.3f, p=%.5f\n', ...
    stats_int_opp.df, stats_int_opp.tstat, p_int_opp);
