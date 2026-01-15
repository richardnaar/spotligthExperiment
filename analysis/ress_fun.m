function [hz, evecs, evals, covAt, maxcomp, maps, snrR, snrE, snrRs, snrEs, elecx, ressx, ress_ts] = ress_fun(data, EEG, nElec, ozLoc, elocs, tidx, peakfreq, compnr, electrode1, plotresults, startfrom) % 
%% ========================================================================
% RESS Function – Rhythmic Entrainment Source Separation
% ========================================================================
% This function is based on the method introduced in:
%
%   Cohen, M. X., & Gulbinaite, R. (2017).
%   *Rhythmic entrainment source separation: Optimizing analyses of neural
%   responses to rhythmic sensory stimulation.* NeuroImage, 147, 43–56.
%   https://doi.org/10.1016/j.neuroimage.2016.11.036
%
% We gratefully acknowledge the original authors, Mike X Cohen and 
% Rasa Gulbinaite, for making their RESS method and code openly available. 
% This function adapts and extends their approach for our analysis of 
% high-frequency SSVEPs in an affective attention paradigm.
%
% Correspondence regarding the original RESS method:
%   mikexcohen@gmail.com or rasa.gulbinaite@gmail.com
%
% NOTE:
% - This implementation includes code derived from the original RESS 
%   demonstration scripts and incorporates minor modifications to fit 
%   our specific preprocessing and analysis pipeline.
% - Additional helper functions are required, including:
%       - filterFGx.m
%       - EEGLAB toolbox (for plotting scalp topographies)
%
% Please cite the above paper if you use or adapt this function.
% =========================================================================


if ~nElec; nElec = EEG.nbchan; end
if ~ozLoc; ozLoc = strcmpi(electrode1,{elocs.labels}); end % ozLoc = strcmpi(electrode1,{EEG.chanlocs.labels})
if isempty(elocs); elocs = EEG.chanlocs; end
ncomps2plot = 4; % this is just for plotting
% set resolution and window parameters
freso =  0.5;
nfft  = ceil( EEG.srate/freso ); % .1 Hz resolution
hz    = linspace(0,EEG.srate,nfft);

peakwidt  = .5; % FWHM at peak frequency
neighfreq = 1;  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;  

% compute covariance matrix at peak frequency 1
fdatAt = filterFGx(data,EEG.srate,peakfreq,peakwidt);
fdatAt = reshape( fdatAt(:,tidx(1):tidx(2),:), nElec,[] );
fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
covAt  = (fdatAt*fdatAt')/diff(tidx);


% compute covariance matrix for lower neighbor 1
fdatLo = filterFGx(data,EEG.srate,peakfreq+neighfreq,neighwidt);
fdatLo = reshape( fdatLo(:,tidx(1):tidx(2),:), nElec,[] );
fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
covLo  = (fdatLo*fdatLo')/diff(tidx);

% compute covariance matrix for upper neighbor
fdatHi = filterFGx(data,EEG.srate,peakfreq-neighfreq,neighwidt);
fdatHi = reshape( fdatHi(:,tidx(1):tidx(2),:), nElec,[] );
fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
covHi  = (fdatHi*fdatHi')/diff(tidx);

% perform generalized eigendecomposition. This is the meat & potatos of RESS
[evecs,evals] = eig(covAt,(covHi+covLo)/2);
[evals,sidx]  = sort(diag(evals),'descend' );
evecs = evecs(:,sidx);

maxcomp = (1:ncomps2plot);

evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)

% maps for plotting
maps = covAt * evecs / (evecs' * covAt * evecs);
[~,idx] = max(abs(maps(:,maxcomp(compnr)))); % find biggest component
maps = maps * sign(maps(idx,maxcomp(compnr))); % force to positive sign        
 

% reconstruct RESS component time series
 ress_ts = zeros(size(data, 2),size(data, 3));
 for ti=1:size(data,3)
     if length(size(data)) < 3
         ress_ts = zeros(size(data,2),1);
         ress_ts(:) = evecs(:,maxcomp(compnr))'*data; % 
     else  
        ress_ts(:,ti) = evecs(:,maxcomp(compnr))'*squeeze(data(:,:,ti)); % 
     end     
 end  
 
% for SNR
ressx = mean(abs( fft(ress_ts(tidx(1):tidx(2),:),nfft,1)/diff(tidx) ).^2,2);
% for elecs
dataX = mean(abs( fft(data(:,tidx(1):tidx(2),:),nfft,2)/diff(tidx) ).^2,3);
elecx = dataX(ozLoc,:,:);


% SNR
[snrR,snrE] = deal(zeros(size(hz)));
skipbins =  10; %  hard-coded!
numbins  = 20+skipbins; %   also hard-coded!

skipbins =  2; %  hard-coded!
numbins  = 4+skipbins; %   also hard-coded!

% loop over frequencies and compute SNR
for hzi=numbins+1:length(hz)-numbins-1
        numer = ressx(hzi);
        denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snrR(hzi) = numer./denom;

        numer = elecx(hzi);
        denom = mean( elecx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snrE(hzi) = numer./denom;
end

numbOfharmonics = 4; fbase = [1:numbOfharmonics];
% choose freq debendent peaks for SNR
hzindex= dsearchn(hz', peakfreq*fbase'); % freqi/2 kui harm

% take the mean around the peak 

snrRs = snrR(hzindex);
snrEs = snrE(hzindex);

if plotresults == 1

    if startfrom == 0
    figure(1), clf
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1); 
    end
    xlim = [3 30];

    cter = startfrom + 1; subplot(2,7,cter)
    plot(hz,snrR,'ro-','linew',1,'markersize',5,'markerface','w')
    hold on
    plot(hz,snrE,'ko-','linew',1,'markersize',5,'markerface','w')
    set(gca,'xlim',xlim)
    axis square
    xlabel('Frequency (Hz)'), ylabel('SNR (1st comp)')

    cter = cter + 1; subplot(2,7,cter)
    map2plot = maps(:,maxcomp(1));
    topoplot(map2plot./max(map2plot),elocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
    title([ 'RESS for (comp 1) ' num2str(peakfreq) ' Hz' ])

    cter = cter + 1; subplot(2,7,cter)
    map2plot = maps(:,maxcomp(2));
    topoplot(map2plot./max(map2plot),elocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
    title([ 'Comp 2 of ' num2str(peakfreq) ' Hz' ])

    cter = cter + 1; subplot(2,7,cter)
    map2plot = maps(:,maxcomp(3));
    topoplot(map2plot./max(map2plot),elocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
    title([ 'Comp 3 of ' num2str(peakfreq) ' Hz' ])

    cter = cter + 1; subplot(2,7,cter)
    map2plot = maps(:,maxcomp(4));
    topoplot(map2plot./max(map2plot),elocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off','shading','interp');
    title([ 'Comp 4 of ' num2str(peakfreq) ' Hz' ])

    cter = cter + 1; subplot(2, 7, cter)
    imagesc(covAt)    
    title([ 'Signal cov (' num2str(size(data, 3)) ')' ])

    cter = cter + 1; subplot(2, 7, cter)
    imagesc(evals)    
    title('Eig vals')

end
end