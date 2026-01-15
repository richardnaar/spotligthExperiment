% Adamian 2019

clear all

mainPath    =   'D:\Adamian2019\';

%locfile     =   'C:\Users\pcadmin\Documents\dok\TARU\eeglab2023.0\64_4EOG.ced';                                     % dir of the electrode location file 
locfile     =   'D:\Adamian2019\64_4EOG_Adamian.ced';
output = 'D:\Adamian2019\output\';

eegPath = [mainPath, 'EEG\'];
eegPathEpoched = [mainPath, 'Epoched data\'];

addpath(mainPath)

eeglab(), close

%% Import


impdir = eegPath; implist = dir([impdir, '*.bdf']);

grandSNR = []; grandAmp = []; grandTrials = [];

subi = 1;
for subi = 1:length(implist)


subName = implist(subi).name;

fprintf('loading participant: %s \n', subName);                                                                     % print the file name to the Command Window

ALLEEG = []; EEG = []; CURRENTSET = [];    

EEG = pop_biosig([impdir, implist(subi).name]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, ...
    'setname', implist(subi).name);   

%% EOG electrodes
exg = {'EXG1'}; 

exgi = find(contains({EEG.chanlocs.labels}, exg));

EEG2 = EEG;
EEG2 = pop_eegfiltnew(EEG2, 1, 0);
EOG = EEG2.data(exgi:exgi+3,:);


EEG.data = EEG.data(1:exgi-1, :); EEG.chanlocs = EEG.chanlocs(1:exgi-1); EEG.nbchan = size(EEG.data,1);

% Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal. 
% EEG = clean_flatlines(EEG, 5);

% EEG = pop_eegfiltnew(EEG, 1, 0);   % high-pass at 1 Hz
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% see EEG.chaninfo.removedchans 
[EEG, rej] = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',8,'norm','on','measure','kurt'); 
    
%% compute average reference
% eegplot(EEG.data, 'events', EEG.event) 
EEG_ref = mean(EEG.data(1:end,:),1); % 1:64

% not adding EOG channel back at the moment
% EEG.data = [EEG.data; EOG]; EEG.chanlocs(end+1) = EOGlocs; EEG.nbchan = EEG.nbchan + 1;

EEG.data = bsxfun(@minus, EEG.data(1:end,:), squeeze(EEG_ref)); % average reference (excluding EOG from the mean)


fprintf(['Averge reference (excluding EOG from the mean) \n'])

%% EOG stuff

% vEOG2 = EOG(1,:)-EOG(2,:);
% hEOG2 = EOG(3,:)-EOG(4,:);

% vEOG2_vel = []; hEOG2_vel = []; 
% for ind=1:length(vEOG2)-1
%     vEOG2_vel(ind) = vEOG2(ind)/vEOG2(ind+1);
%     hEOG2_vel(ind) = hEOG2(ind)/hEOG2(ind+1);
% end

%
% seg_dur = 3*EEG.srate;
% transient = floor(0.4*EEG.srate);
% count_trials = 0; count_out_trials = 0;
% for eventi = 1:length(EEG.event)
%     if any(EEG.event(eventi).type == 1:6)
%         count_trials = count_trials + 1;
%         grandTrials(subi) = count_trials;
%         current_data = vEOG2(:, EEG.event(eventi).latency+transient:EEG.event(eventi).latency+seg_dur);
%         current_data = current_data-mean(current_data(1:transient));
%         maxv = (max(abs(current_data)));
% 
%         current_data2 = hEOG2(:, EEG.event(eventi).latency+transient:EEG.event(eventi).latency+seg_dur);
%         current_data2 = current_data2-(mean(current_data2(1:transient)));
% 
%         maxh = (max(abs(current_data2)));
%         if maxh > median(current_data2)+5*std(current_data2)
%             count_out_trials = count_out_trials + 1;
%             % plot(current_data2')
%             % pause()
%         end
%         plot(current_data')
%         hold on
%         plot(current_data2')
%         title(['maxv: ' num2str(maxv) ',  maxh: ' num2str(maxh)]);
%         legend('vEOG','hEOG')
%         pause()
%         hold off
%     end
% end

%% triggers

% Trigger codes 1-6 mark the start of each trial with the following condition mapping:
% Trigger code 1: attend left blue --cyan
% Trigger code 2: attend reght red -pink
% Trigger code 3: attend left red --pink
% Trigger code 4: attend right blue -cyan
% Trigger code 5: attend both blue -blue
 


for condi = 1:6
    counter = 1;
    onset = []; 
    for eventi = 1:length(EEG.event)
        if EEG.event(eventi).type == condi % < 7 
            onset(counter) = EEG.event(eventi).latency;
            counter = counter + 1;
        end
    end

%% 
pool = {'PO3','PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'O1', 'O2', 'I1', 'I2', 'Oz', 'Iz', 'POz'};
% pool = {'PO3'};
% electrode = 'PO3'; % 'PO8';
elec2plot = find(ismember({EEG.chanlocs.labels}, pool));


startInS = 0.4;
dur = 2.9;
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

% hz = 0:1/(dur-startInS):EEG.srate/2;
% hz = (0:length(juku)-1) * (EEG.srate / length(juku));


% stem(hz(2:500), juku(2:500))
% plot(hz(2:100), juku(2:100))

% eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

skipbins = 1; 
numbins = 1 + skipbins;

% Computing SNR
for hzi = numbins + 1:length(juku) - numbins - 1
    numer = juku(hzi);
    denom = rms(juku([hzi-numbins:hzi-skipbins, hzi+skipbins:hzi+numbins]));
    snrE(hzi) = numer./denom;
end 

% plot(hz(2:100), snrE(2:100))

grandAmp(subi).cond(condi).amp = juku; % save amplitude instead
grandSNR(subi).cond(condi).snr = snrE; 

end
end

save([output, 'grandSNR-', date, '.mat'], 'grandSNR', '-v7.3');
save([output, 'grandAmp-', date, '.mat'], 'grandAmp', '-v7.3');

%% prepare data for plotting
output_list = dir([output, '*.mat']);
load([output, output_list(1).name])

hz = 0:1/(dur-startInS):EEG.srate/2;

% meanSNR = mean(cat(1, grandSNR.snr), 1);
current_data = grandAmp; % grandSNR; % or grandAmp

nSubj = numel(current_data);
nCond = numel(current_data(1).cond);

% Find how many timepoints
fbin = numel(current_data(1).cond(1).amp); % or amp

% Preallocate [frex × cond × subj]
datForPlotting = nan(fbin, nCond, nSubj);

for s = 1:nSubj
    for c = 1:nCond
        datForPlotting(:, c, s) = current_data(s).cond(c).amp(:); % or amp
    end
end

avg_dat = mean(datForPlotting,3);
hz_inds = dsearchn(hz', [7.5,10,12, 8.57]');
hz(hz_inds)
avg_dat = avg_dat(hz_inds,:);

%% Plot the raw SNR
% Trigger code 1: attend left blue --cyan
% Trigger code 2: attend reght red -magenta
% Trigger code 3: attend left red --magenta
% Trigger code 4: attend right blue -cyan
% Trigger code 5: attend both blue -blue
% Trigger code 6: attend both red -red

cols = {'--cyan','-magenta','--magenta','-cyan','blue','red'};
conOrder = [1,4,5,3,2,6];

figure(1); hold
for plti = 1:length(cols)
    current_data = mean(datForPlotting(:, conOrder(plti), :),3, 'omitnan');
    plot(hz(6:50), current_data(6:50), cols{conOrder(plti)})
end

ylabel('Amplitude', 'FontSize',14); xlabel('Frequency (Hz)', 'FontSize',14)  

set(findall(gcf,'-property','FontSize'),'FontSize',14)

legend({'attend left blue';'attend right blue';'attend both blue';'attend left red';'attend right red';'attend both red'})
% legend({'attend left blue';'attend right red';'attend left red';'attend right blue';'attend both blue';'attend both red'})

%% Bar plot

% Define colors: first 3 blue, last 3 red
colors = [repmat([0 0 1], 3, 1); repmat([1 0 0], 3, 1)];

N = 16;
err = avg_dat./sqrt(N);
conOrder2 = [1,5,4,3,6,2];

% Plot
figure;
for f = 1:size(avg_dat, 1)
    % subplot(2, 2, f);  % adjust layout as needed
    data = avg_dat(f, conOrder2);
    e = err(f, conOrder2);

    b = bar(data, 'FaceColor', 'flat');
    for c = 1:length(conOrder2)
        b.CData(c,:) = colors(c,:);
    end
    
    % title(sprintf('Frequency %d', f));
    ylabel('Amplitude');
    xlabel('Condition');
    xticks(1:6);
    xticklabels(conOrder2);
    ylim([0 max(avg_dat(:))*1.4]);
    
    hold on
    errorb(1:6,data,e,'top')

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','fontname'),'fontname', 'Arial')
    pause()
    hold off
end

%% Epoched (turns out that there is no trial rejection information)

cleanlist = dir([eegPathEpoched,'*.set']);     

subi = 2;

% load raw EEG dataset
ALLEEG = []; EEG = []; CURRENTSET = 1;
EEG = pop_loadset('filename', cleanlist(subi).name, 'filepath', eegPathEpoched);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

% load raw EOG dataset
EOG = pop_loadset('filename', cleanlist(subi+1).name, 'filepath', eegPathEpoched);

%% ADD COMPLEX DEMOD

%% Parameters
params.targetFreqs   = [7.5, 10, 12, 8.57];   % Hz; SSVEP frequencies of interest
params.startInS      = 0.4;      % seconds after trigger to start analysis window
params.durInS        = 2.9;      % length of analysis window in seconds

% Taper (square–cosine ramp) as fraction of window length at each edge
% e.g. 0.1 = 10% of the samples at each edge tapered
params.taperFraction = 0.1;      

% Low-pass for envelope (how “slowly” amplitude can vary)
params.lpCutoffHz    = 1.0;      % Hz (0.5–2 Hz is typical for SSVEP envelope)
params.lpOrder       = 3;        % Butterworth order (3–6 is usually fine)

% Optional baseline (if you later want baseline correction of demod amplitude)
params.baseline      = [];       % e.g. [-0.5 0];  % relative to trigger, in seconds
params.normflag      = 0;        % 0=subtractive, 1=divisive baseline (if used)

params.epochPre   = -0.5;   % seconds relative to event time (e.g. sync or flicker onset)
params.epochPost  =  3.5;   % seconds after event

% analysis window relative to *flicker onset* (adapt to paper)
% If your event is at flicker onset, you can use:
params.analysisStart = 0.526;  % s after flicker onset
params.analysisEnd   = 2.9;% 3.026;  % s after flicker onset
% only for RESS
params.compnr        = 1;    
params.ressSegmFilter = 10;
%% Import


impdir = eegPath; implist = dir([impdir, '*.bdf']);

grandSNR = []; grandAmp = []; grandTrials = [];
grandSNR_ress = []; grandAmp_ress = [];
% subi = 3;
for subi = 1:length(implist)


subName = implist(subi).name;

fprintf('loading participant: %s \n', subName);                                                                     % print the file name to the Command Window

ALLEEG = []; EEG = []; CURRENTSET = [];    

EEG = pop_biosig([impdir, implist(subi).name]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, ...
    'setname', implist(subi).name);   

%% Cant use this for Adamian2019
% EEG = pop_chanedit(EEG, 'load',{locfile 'filetype' 'autodetect'});                                     % edit channel locations structure of an EEGLAB dataset, EEG.chanlocs. 
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

%% For Adamian2019 use instead
locs = readlocs(locfile);

% Force assignment by channel index
for ch = 1:min(length(EEG.chanlocs), length(locs))
    EEG.chanlocs(ch).X          = locs(ch).X;
    EEG.chanlocs(ch).Y          = locs(ch).Y;
    EEG.chanlocs(ch).Z          = locs(ch).Z;
    EEG.chanlocs(ch).theta      = locs(ch).theta;
    EEG.chanlocs(ch).radius     = locs(ch).radius;
    EEG.chanlocs(ch).sph_theta = locs(ch).sph_theta;
    EEG.chanlocs(ch).sph_phi   = locs(ch).sph_phi;
    EEG.chanlocs(ch).sph_radius= locs(ch).sph_radius;
end

EEG = eeg_checkset(EEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
%%
%% EOG electrodes
exg = {'EXG1'}; 

exgi = find(contains({EEG.chanlocs.labels}, exg));

EEG2 = EEG;
EEG2 = pop_eegfiltnew(EEG2, 1, 0);
EOG = EEG2.data(exgi:exgi+3,:);


EEG.data = EEG.data(1:exgi-1, :); EEG.chanlocs = EEG.chanlocs(1:exgi-1); EEG.nbchan = size(EEG.data,1);

% Maximum tolerated flatline duration. In seconds. If a channel has a longer flatline than this, it will be considered abnormal. 
% EEG = clean_flatlines(EEG, 5);

% EEG = pop_eegfiltnew(EEG, 1, 0);   % high-pass at 1 Hz
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 

% see EEG.chaninfo.removedchans 
[EEG, rej] = pop_rejchan(EEG, 'elec', 1:size(EEG.data,1) ,'threshold',8,'norm','on','measure','kurt'); 
    
%% compute average reference
% eegplot(EEG.data, 'events', EEG.event) 
EEG_ref = mean(EEG.data(1:end,:),1); % 1:64

% not adding EOG channel back at the moment
% EEG.data = [EEG.data; EOG]; EEG.chanlocs(end+1) = EOGlocs; EEG.nbchan = EEG.nbchan + 1;

EEG.data = bsxfun(@minus, EEG.data(1:end,:), squeeze(EEG_ref)); % average reference (excluding EOG from the mean)


fprintf(['Averge reference (excluding EOG from the mean) \n'])

%% EOG stuff

% vEOG2 = EOG(1,:)-EOG(2,:);
% hEOG2 = EOG(3,:)-EOG(4,:);

% vEOG2_vel = []; hEOG2_vel = []; 
% for ind=1:length(vEOG2)-1
%     vEOG2_vel(ind) = vEOG2(ind)/vEOG2(ind+1);
%     hEOG2_vel(ind) = hEOG2(ind)/hEOG2(ind+1);
% end

%
% seg_dur = 3*EEG.srate;
% transient = floor(0.4*EEG.srate);
% count_trials = 0; count_out_trials = 0;
% for eventi = 1:length(EEG.event)
%     if any(EEG.event(eventi).type == 1:6)
%         count_trials = count_trials + 1;
%         grandTrials(subi) = count_trials;
%         current_data = vEOG2(:, EEG.event(eventi).latency+transient:EEG.event(eventi).latency+seg_dur);
%         current_data = current_data-mean(current_data(1:transient));
%         maxv = (max(abs(current_data)));
% 
%         current_data2 = hEOG2(:, EEG.event(eventi).latency+transient:EEG.event(eventi).latency+seg_dur);
%         current_data2 = current_data2-(mean(current_data2(1:transient)));
% 
%         maxh = (max(abs(current_data2)));
%         if maxh > median(current_data2)+5*std(current_data2)
%             count_out_trials = count_out_trials + 1;
%             % plot(current_data2')
%             % pause()
%         end
%         plot(current_data')
%         hold on
%         plot(current_data2')
%         title(['maxv: ' num2str(maxv) ',  maxh: ' num2str(maxh)]);
%         legend('vEOG','hEOG')
%         pause()
%         hold off
%     end
% end

%% triggers

% Trigger codes 1-6 mark the start of each trial with the following condition mapping:
% Trigger code 1: attend left blue --cyan
% Trigger code 2: attend reght red -pink
% Trigger code 3: attend left red --pink
% Trigger code 4: attend right blue -cyan
% Trigger code 5: attend both blue -blue
 
pool = {'PO3','PO4', 'PO7', 'PO8', 'PO9', 'PO10', 'O1', 'O2', 'I1', 'I2', 'Oz', 'Iz', 'POz'};

% ress_pool = { ...
% 'O1','O2','Oz','Iz', 'I1','I2', ...
% 'PO3','PO4','PO7','PO8','PO9','PO10','POz', ...
% 'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','Pz', ...
% 'CP1','CP2','CP3','CP4','CP5','CP6','CPz', ...
% 'TP7','TP8' };

% pool = {'PO3'};
% electrode = 'PO3'; % 'PO8';
elec2plot = find(ismember({EEG.chanlocs.labels}, pool));
elec_ress = find(ismember({EEG.chanlocs.labels}, pool));
%% Compute RESS
plotresults = 0; electrode = 'Oz'; startfrom = 0; % just for plotting, not using

% Set up a data matrix
srate = EEG.srate;
data = []; % zeros(size(EEG.chanlocs,2), round(srate*(params.durInS+2)), length(EEG.event)/2);
cond_info_ress_ts = [];
trialCounter = 0;
for eind = 1:length(EEG.event)
    if any(EEG.event(eind).type == [1,2,3,4,5,6])
        trialCounter = trialCounter + 1;
        % 1 second befor and 1 after the end of the event
        try
            % data(:,:,trialCounter) =  EEG.data(:, floor(EEG.event(eind).latency-srate):floor(EEG.event(eind).latency+srate*(params.durInS+1))-1);
            % data(:,:,trialCounter) =  EEG.data(elec_ress, floor(EEG.event(eind).latency-srate):floor(EEG.event(eind).latency+srate*(params.analysisEnd+1))-1);
            data(:,:,trialCounter) =  EEG.data(elec_ress, floor(EEG.event(eind).latency-srate):floor(EEG.event(eind).latency+srate*(params.ressSegmFilter-1))-1); % 10 s window for filtering / hard coded!

            cond_info_ress_ts(trialCounter) = EEG.event(eind).type; 
        catch
            if eind == length(EEG.event)
                fprintf('Error with segmentation: have to leave the last trial out!\n')
            else
                fprintf('Unidentified segmentation error\n')
            end
        end
    end
end

% train ress with data taken from 0.6 s after the event (1.6 s after epoch start)
% tidx = [round(srate*1.6), size(data,2)-srate]; % we use data from start + 0.6 to -6.1 seconds (epoch starts 1s before the event and ends 1s after the end of the stimulation)
tidx = [round(srate*1.6), size(data,2)-round(srate*(params.ressSegmFilter-3.9))]; % we use data from start + 0.6 to -6.1 seconds (epoch starts 1s before the event and ends 1s after the end of the stimulation)
ress_ts = zeros(length(params.targetFreqs), size(data,2),size(data,3));
for freqi = 1:length(params.targetFreqs)
    % output: hz, evecs, evals, covAt, maxcomp, maps, snrR, snrE, snrRs, snrE, elecx, ressx, ress_time_series
    [~, evecs, ~,  ~, maxcomp, maps, ~, ~, ~, ~, ~, ~, ~] = ...
              ress_fun(mean(data,3), EEG, length(elec_ress), 0, EEG.chanlocs(elec_ress), tidx, params.targetFreqs(freqi), params.compnr, electrode, plotresults, startfrom);

    ress_time_series = evecs2ress(data, evecs, 1); % time x trials
    ress_ts(freqi,:,:) = ress_time_series;
end
%%
% ress_ts = squeeze(sum(ress_ts,1)); % this is probably not helpful

for condi = 1:6

onset = [EEG.event([EEG.event.type] == condi).latency]; 

meanComplexFFT = [];
for eventi = 1:length(onset)
    transient = params.startInS *EEG.srate; % in samples
    eventDur = params.durInS*EEG.srate;
    % current_data = EEG.data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) );
    current_data = mean(EEG.data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) ),1); % avg over electrodes
    
    fft_dat = fft(current_data);
    
    meanComplexFFT(:,eventi) = fft_dat/length(current_data); % coherent

    % meanComplexFFT(:,eventi) = abs(fft_dat)/length(current_data)*2;
end

juku = abs(mean(meanComplexFFT,2)); % lose coherence
% juku = mean(meanComplexFFT,2);

% hz = 0:1/(size(meanComplexFFT,1)/EEG.srate):EEG.srate/2;
% hz = (0:length(juku)-1) * (EEG.srate / length(juku));


% stem(hz(2:500), juku(2:500))
% plot(hz(2:100), juku(2:100))

% eegplot(EEG.data, 'srate', EEG.srate, 'events', EEG.event);

snrE = nan(size(juku));
skipbins = 1; 
numbins = 1 + skipbins;

% Computing SNR
for hzi = numbins + 1:length(juku) - numbins - 1
    numer = juku(hzi);
    denom = rms(juku([hzi-numbins:hzi-skipbins, hzi+skipbins:hzi+numbins]));
    snrE(hzi) = numer./denom;
end 

% plot(hz(2:100), snrE(2:100))

grandAmp(subi).cond(condi).amp = juku; % save amplitude instead
grandSNR(subi).cond(condi).snr = snrE; 

%% RESS

cond_idx = find(cond_info_ress_ts == condi);
transient = 1+params.startInS*EEG.srate; % in samples, the RESS epoch starts 1s before
eventDur = params.durInS*EEG.srate;


for f = 1:size(ress_ts,1)

    meanComplexRESS = [];
    for eventi = 1:length(cond_idx)
         tr = cond_idx(eventi);
        % current_data = mean(data(elec2plot, round(onset(eventi)+transient):round((onset(eventi)+(eventDur)-1)) ),1); % avg over electrodes
        % current_data = ress_ts(f, round(transient):round(eventDur), eventi);
        current_data = ress_ts(f, round(transient):round(eventDur), tr);
    
        fft_dat = fft(current_data);
        
        meanComplexRESS(:,eventi) = fft_dat/length(current_data); % coherent
    
        % meanComplexFFT(:,eventi) = abs(fft_dat)/length(current_data)*2;
    end
    
    ress_amp = abs(mean(meanComplexRESS,2)); % lose coherence
    
    snrRESS = nan(size(ress_amp));  % <-- reset!
    skipbins = 1; 
    numbins = 1 + skipbins;
    
    % Computing SNR
    for hzi = numbins + 1:length(ress_amp) - numbins - 1
        numer = ress_amp(hzi);
        denom = rms(ress_amp([hzi-numbins:hzi-skipbins, hzi+skipbins:hzi+numbins]));
        snrRESS(hzi) = numer./denom;
    end 
    
    hz = 0:1/(size(ress_amp,1)/EEG.srate):EEG.srate/2;
    % plot(hz(2:100), snrRESS(2:100))
    % plot(hz(2:100), ress_amp(2:100))
    
    grandAmp_ress(subi).cond(condi).amp = ress_amp; % save amplitude instead
    grandSNR_ress(subi).cond(condi).snr = snrRESS; 

end

%% complex demod 
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
end



save([output, 'grandSNR-', date, '.mat'], 'grandSNR', '-v7.3');
save([output, 'grandAmp-', date, '.mat'], 'grandAmp', '-v7.3');
save([output, 'grandSNR_ress-', date, '.mat'], 'grandSNR_ress', '-v7.3');
save([output, 'grandAmp_ress-', date, '.mat'], 'grandAmp_ress', '-v7.3');
save([output, 'grandDemod-', date, '.mat'], 'grandDemod', '-v7.3');

%% Grand Avg across frex x cond

% % add loading here
% data2plot = grandAmp_ress;
% 
% %% prepare data for plotting
% output_list = dir([output, '*.mat']);
% load([output, output_list(1).name])
% 
% startInS = 0.4;
% dur = 2.9;
% 
% hz = 0:1/(dur-startInS):EEG.srate/2;
% 
% % meanSNR = mean(cat(1, grandSNR.snr), 1);
% current_data = grandAmp; % grandSNR; % or grandAmp
% 
% nSubj = numel(current_data);
% nCond = numel(current_data(1).cond);
% 
% % Find how many timepoints
% fbin = numel(current_data(1).cond(1).amp); % or amp
% 
% % Preallocate [frex × cond × subj]
% datForPlotting = nan(fbin, nCond, nSubj);
% 
% for s = 1:nSubj
%     for c = 1:nCond
%         datForPlotting(:, c, s) = current_data(s).cond(c).amp(:); % or amp
%     end
% end
% 
% avg_dat = mean(datForPlotting,3);
% hz_inds = dsearchn(hz', [7.5,10,12, 8.57]');
% hz(hz_inds)
% avg_dat = avg_dat(hz_inds,:);

%% 
outlist = dir([output, '*.mat']);

load([output, outlist(4).name]) % 4,8


%% Skip if not demod

targetFreqs = [7.5 10 12 8.57];

srate = 256;

grandSNR_ress_scalar = struct([]);

data_current = grandSNR; % grandAmp;% grandAmp_ress; % grandSNR_ress

nSub  = numel(data_current);
nCond = numel(data_current(1).cond);

for subi = 1:nSub
    for condi = 1:nCond

        % snr = data_current(subi).cond(condi).snr(:);  % N x 1
        snr = data_current(subi).cond(condi).amp(:);  % N x 1
        N   = numel(snr);

        % Full-spectrum Hz vector matching snr length exactly:
        hz_full = (0:N-1) * (srate / N);

        % Indices of bins closest to your target frequencies:
        hz_inds = dsearchn(hz_full', targetFreqs');

        % Extract 4 scalar values (like complex demod meanAmp_scalar):
        grandSNR_ress_scalar(subi).cond(condi).meanAmp_scalar = snr(hz_inds);

        % Store frequency list (your plotting code uses this)
        grandSNR_ress_scalar(subi).cond(condi).freqs = targetFreqs;

        % (Optional) store for debugging / sanity checks
        grandSNR_ress_scalar(subi).cond(condi).hz_inds = hz_inds;
        grandSNR_ress_scalar(subi).cond(condi).hz_at   = hz_full(hz_inds);
    end
end

% Now you can do:
data2plot = grandSNR_ress_scalar;


%% Skip end

%NB!
data2plot = grandDemod; % grandSNR_ress;% grandDemod;


% ----- collect amplitudes from grandDemod into an array -----
nSub  = numel(data2plot);
nCond = numel(data2plot(1).cond);
nFreq = numel(data2plot(1).cond(1).meanAmp_scalar);

allAmp = nan(nSub, nCond, nFreq);   % subj × cond × freq

for subi = 1:nSub
    for condi = 1:nCond
        allAmp(subi,condi,:) = data2plot(subi).cond(condi).meanAmp_scalar(:);
    end
end

% ----- mean and SEM across subjects -----
avg_dat = squeeze(mean(allAmp,1));          % cond × freq
avg_dat = avg_dat.';                        % freq × cond (to match your old code)

N   = nSub;
err = squeeze(std(allAmp,[],1)).'/sqrt(N);  % freq × cond  (SEM)

freqs = data2plot(1).cond(1).freqs; 

%%

% Define colors: first 3 blue, last 3 red
colors = [repmat([0 0 1], 3, 1); repmat([1 0 0], 3, 1)];

conOrder2 = [1,5,4,3,6,2];   % your preferred condition order

figure;
for f = 1:size(avg_dat, 1)    % loop over frequencies

    data = avg_dat(f, conOrder2);
    e    = err(f,  conOrder2);

    b = bar(data, 'FaceColor', 'flat');
    hold on
    for c = 1:length(conOrder2)
        b.CData(c,:) = colors(c,:);
    end
    
    ylabel('Complex demod');
    xlabel('Condition');
    xticks(1:6);
    xticklabels(conOrder2);
    ylim([0 max(avg_dat(:))*1.4]);
    
    % error bars (using your errorb function)
    errorb(1:6, data, e, 'top');

    % optional: add frequency in the title
    % freqs(f) assumed from above
    title(sprintf('Frequency %.1f Hz', freqs(f)));

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(findall(gcf,'-property','fontname'),'fontname', 'Arial')

    pause();   % step through frequencies
    hold off
end


%% Delete

nSub  = length(grandDemod);
nCond = 6;
nFreq = length(grandDemod(1).cond(1).freqs);

allAmp = nan(nSub, nCond, nFreq);

for subi = 1:nSub
    for condi = 1:nCond
        allAmp(subi, condi, :) = grandDemod(subi).cond(condi).meanAmp_scalar;
    end
end

% mean across subjects
grandAmp_byCond = squeeze(nanmean(allAmp, 1));  
% size: nCond x nFreq

%% Delete

nTime = length(grandDemod(1).cond(1).tEpoch);

allAmp_time = nan(nSub, nCond, nFreq, nTime);

for subi = 1:nSub
    for condi = 1:nCond
        allAmp_time(subi, condi, :, :) = grandDemod(subi).cond(condi).demAmpFreq;
    end
end

% average across subjects
grandAmpTime_byCond = squeeze(nanmean(allAmp_time, 1));

condi = 3;
fi = 1;   % first frequency

plot(grandDemod(1).cond(1).tEpoch, squeeze(grandAmpTime_byCond(condi,fi,:)));
xlabel('Time (s)');
ylabel('SSVEP amplitude');
title(sprintf('Condition %d, Frequency %.1f Hz', condi, targetFreqs(fi)));

%% cronbach
cronbach_cond = [];
for f = 1:4
    cronbach_cond(f) = cronbach(squeeze(allAmp(:,:,f)));
end

mean(cronbach_cond) 
% 
% demod:     0.9547 (0.9577    0.9463    0.9353    0.9797) 
% coherence: 0.8997 (0.8713    0.9364    0.9555    0.8354)
% ress:      0.8796 (0.8354    0.8251    0.9223    0.9358)
%%
% collapse across conditions first, then average across subjects
amp_collapsedCond = squeeze(nanmean(allAmp, 2));
% size: nSub x nFreq

grandAmp_overall = nanmean(amp_collapsedCond, 1);
% size: 1 x nFreq

%% table

Subject = [];
Condition = [];
Frequency = [];
Amplitude = [];

for subi = 1:nSub
  for condi = 1:nCond
    for fi = 1:nFreq
        Subject(end+1)   = subi;
        Condition(end+1) = condi;
        Frequency(end+1) = targetFreqs(fi);
        Amplitude(end+1) = grandDemod(subi).cond(condi).meanAmp_scalar(fi);
    end
  end
end

grandTable = table(Subject', Condition', Frequency', Amplitude', ...
                   'VariableNames', {'Sub','Cond','Freq','Amp'});