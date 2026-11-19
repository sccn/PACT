% pac_scanHfoByPhase() - It generates a matrix of phase vs. amp which
%                        allows one to search for parameter space.
%
% Usage:
%        EEG = pac_scanHfoByPhase(EEG, phaseFreqRange, ampFreqRange, ...
%                                 HASnumber, numFreqBins, plotType, ...
%                                 normalizeColor)
% Input:
%        EEG -- EEGLAB data structure.
%        phaseFreqRange -- [lowHz highHz]
%        ampFreqRange -- [lowHz highHz]
%        HASnumber -- [%]
%        numFreqBins -- [integer]
%        plotType -- [1|2|3] where 3 disables plotting (i.e. batch mode)
%        normalizeColor -- [0|1] color scale normalization.
%
%
%

% Author: Makoto Miyakoshi JSPS/SCCN,INC,UCSD
% History:
% 01/13/2020 Makoto. Check before moving to Github.
% 11/29/2020 Makoto. Rayleigh test added.
% 11/22/2020 Makoto. Fixed bug related to hilbert() input--must be column-wise operation.
% 11/20/2020 Makoto. Debugged on analyticPhase calculation.
% 11/07/2020 Makoto. Implemented Daisuke's idea.
% 10/30/2020 Makoto. Updated. Batch mode supported.
% 07/24/2019 Makoto. Created.

% Copyright (C) 2019, Makoto Miyakoshi SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function EEG = pac_scanHfoByPhase(EEG, phaseFreqRange, ampFreqRange, HASnumber, numFreqBins, plotType, normalizeColor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add the necessary folders to the Matlab path. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('circ_rtest','file')
    p = which('eegplugin_pac');
    p = p(1:findstr(p,'eegplugin_pac.m')-1);
    p = genpath(p);
    addpath(p);
end

if EEG.trials ~= 1 || ndims(EEG.data) ~= 2
    error('PACT scanning supports continuous data only.');
end
if numFreqBins < 4 || numFreqBins ~= round(numFreqBins) || HASnumber <= 0 || HASnumber > 100 || ...
        any(phaseFreqRange<=0) || any(ampFreqRange<=0) || ...
        any(diff(phaseFreqRange)<=0) || any(diff(ampFreqRange)<=0) || ampFreqRange(2)>=EEG.srate/2
    error('Invalid frequency scan parameters.');
end

if isfield(EEG, 'pacScan')
    EEG = rmfield(EEG,'pacScan');
end

% Prepare log-scaled frequency edges.
phaseFreqList = logspace(log10(phaseFreqRange(1)), log10(phaseFreqRange(2)), numFreqBins);
ampFreqList   = logspace(log10(ampFreqRange(1)),   log10(ampFreqRange(2)),   numFreqBins);
nPhase = numel(phaseFreqList)-1;
nAmp   = numel(ampFreqList)-1;

% Hamming windowed-sinc filters (firfilt). The historical transition widths
% were one value per family (phase: half the mean bin width; amplitude: the
% mean bin width). For log-spaced bins narrower than that width the passband
% collapsed: the centre gain fell below 1 and pop_firws warned "Filter order
% too low". Each band's width is now capped at pop_firws' limit, so the -6 dB
% edges stay at the requested bin edges (see pac_firwsOrder, VALIDATION.md).
phaseFreqTbw = mean(diff(phaseFreqList))/2;
ampFreqTbw   = mean(diff(ampFreqList));
phaseFilterOrder = zeros(1,nPhase);
ampFilterOrder   = zeros(1,nAmp);
for k = 1:nPhase
    phaseFilterOrder(k) = pac_firwsOrder(EEG.srate,phaseFreqList(k:k+1),phaseFreqTbw);
end
for k = 1:nAmp
    ampFilterOrder(k) = pac_firwsOrder(EEG.srate,ampFreqList(k:k+1),ampFreqTbw);
end
fprintf('Phase bands (Hz): %s\n',mat2str(round(phaseFreqList,3)));
fprintf('Phase filter orders: %s\n',mat2str(phaseFilterOrder));
fprintf('Amplitude bands (Hz): %s\n',mat2str(round(ampFreqList,3)));
fprintf('Amplitude filter orders: %s\n',mat2str(ampFilterOrder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exclude filter transients at recording edges and boundaries.    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The guard covers half the longest filter (its transient support) and at
% least 1 s. Hilbert transforms are computed within boundary segments.
workEEG = EEG;
isBoundary = false(1,numel(EEG.event));
for k = 1:numel(EEG.event)
    value = EEG.event(k).type;
    isBoundary(k) = isequal(value,-99) || ...
        ((ischar(value) || (isstring(value) && isscalar(value))) && strcmpi(value,'boundary'));
end
workEEG.event = EEG.event(isBoundary);
for k = 1:numel(workEEG.event)
    workEEG.event(k).type = 'boundary';
end
starts = findboundaries(workEEG.event);
cuts = [starts(starts>=1 & starts<=EEG.pnts) EEG.pnts+1];
guard = max(round(EEG.srate),max([phaseFilterOrder ampFilterOrder])/2);
validSamples = false(1,EEG.pnts);
for s = 1:numel(cuts)-1
    span = cuts(s):cuts(s+1)-1;
    if numel(span) > 2*guard
        validSamples(span(1)+guard:span(end)-guard) = true;
    end
end
boundaryLatency = cuts(2:end-1);
periBoundaryIdx = find(~validSamples);
if ~any(validSamples)
    error('PACT:ScanUnsupported',['No samples remain after excluding %.1f s filter guards ',...
        'at recording edges and boundaries. Use wider or higher frequency bins, or longer ',...
        'continuous data.'],guard/EEG.srate);
end
fprintf('%d boundaries found. Guard %.2f s; %.1f s of data excluded.\n',...
    numel(boundaryLatency),guard/EEG.srate,numel(periBoundaryIdx)/EEG.srate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop phase x amplitude combinations for all the channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HAS threshold is computed separately for each channel and band.
resulVectLen = zeros(EEG.nbchan,nPhase,nAmp);
mAngle       = zeros(EEG.nbchan,nPhase,nAmp);
modIdx       = zeros(EEG.nbchan,nPhase,nAmp);
hfoAmp       = zeros(EEG.nbchan,nPhase,nAmp);
rayleighZ    = zeros(EEG.nbchan,nPhase,nAmp);

% Amplitude envelopes do not depend on the phase band: compute them once.
ampEnvelopes = cell(1,nAmp);
for ampFreqIdx = 1:nAmp
    filtered = pop_firws(workEEG, 'fcutoff', ampFreqList(ampFreqIdx:ampFreqIdx+1), 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', ampFilterOrder(ampFreqIdx), 'minphase', 0);
    ampEnvelopes{ampFreqIdx} = segmentHilbert(filtered.data,cuts,validSamples,'abs');
end
clear filtered

EEGorig = EEG;
for phaseFreqIdx = 1:nPhase

    % Apply band-pass filter to extract low-freq phase.
    filtered = pop_firws(workEEG, 'fcutoff', phaseFreqList(phaseFreqIdx:phaseFreqIdx+1), 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', phaseFilterOrder(phaseFreqIdx), 'minphase', 0);
    analyPhase = segmentHilbert(filtered.data,cuts,validSamples,'angle');

    for ampFreqIdx = 1:nAmp
        analyAmp = ampEnvelopes{ampFreqIdx};
        thresholdVector = prctile(analyAmp', 100-HASnumber);

        for chIdx = 1:EEG.nbchan
            currentChHfoMask = HASnumber==100 | analyAmp(chIdx,:)>thresholdVector(chIdx);
            tmpChHfo   = analyAmp(  chIdx,currentChHfoMask);
            tmpChPhase = analyPhase(chIdx,currentChHfoMask);
            
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %%% Visualize for validation (11/29/2020) %%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if ampFreqIdx ==5 & chIdx == 9
%                 
%                 figure
%                 plot(analyAmp(chIdx,:))
%                 hold on
%                 plot(find(currentChHfoMask), analyAmp(chIdx,find(currentChHfoMask)), 'r.')
%                 
%                 phaseFreqBinCenter = phaseFreqList(1:end-1)+diff(phaseFreqList)/2;
%                 
%                 figure
%                 circ_plot(tmpChPhase','hist',[], 36, true,true,'linewidth',3,'color','r');
%                 suptitle(sprintf('Phase%.2f, RVL %.3f, p=%.3f', phaseFreqBinCenter(phaseFreqIdx), circ_r(tmpChPhase'), circ_rtest(tmpChPhase')))
%             end
            
            if isempty(tmpChPhase)
                rayleighZ(chIdx,phaseFreqIdx,ampFreqIdx) = NaN;
                resulVectLen(chIdx,phaseFreqIdx,ampFreqIdx) = NaN;
                mAngle(chIdx,phaseFreqIdx,ampFreqIdx) = NaN;
                modIdx(chIdx,phaseFreqIdx,ampFreqIdx) = NaN;
                hfoAmp(chIdx,phaseFreqIdx,ampFreqIdx) = NaN;
                continue
            end
            [~, zStats] = circ_rtest(tmpChPhase');
            rayleighZ(   chIdx,phaseFreqIdx,ampFreqIdx) = zStats;
            resulVectLen(chIdx,phaseFreqIdx,ampFreqIdx) = circ_r(tmpChPhase');
            mAngle(      chIdx,phaseFreqIdx,ampFreqIdx) = circ_rad2ang(circ_mean(tmpChPhase'));
            modIdx(      chIdx,phaseFreqIdx,ampFreqIdx) = abs(mean(tmpChHfo.*exp(1i*tmpChPhase)));
            hfoAmp(      chIdx,phaseFreqIdx,ampFreqIdx) = mean(tmpChHfo);
            
            % Per-cell phase samples have unequal lengths. Stacking them
            % by repetition changes the observations; this legacy field
            % has been retired in the scientifically validated output.

        end
    end
end
% close(waitbarHandle)

% Store the results.
EEG = EEGorig;
EEG.pacScan.resulVectLen = single(resulVectLen);
EEG.pacScan.pacAngle     = single(mAngle);
EEG.pacScan.lfoPhaseMi   = single(modIdx);
EEG.pacScan.meanHfoAmp   = single(hfoAmp);
EEG.pacScan.rayleighZ    = single(rayleighZ);
EEG.pacScan.analyticPhaseStack = []; % Deprecated: unequal-length HAS samples cannot be stacked safely.
EEG.pacScan.phaseFreqEdges = phaseFreqList;
EEG.pacScan.ampFreqEdges   = ampFreqList;
EEG.pacScan.phaseFreqTbw   = phaseFreqTbw;
EEG.pacScan.ampFreqTbw     = ampFreqTbw;
EEG.pacScan.phaseFilterOrder = phaseFilterOrder;
EEG.pacScan.ampFilterOrder   = ampFilterOrder;
EEG.pacScan.guardSamples     = guard;
EEG.pacScan.validSampleMask  = validSamples;
EEG.pacScan.dataDimensions = {'Channels', 'Phase Freq Edges [Hz]', 'HFO Freq Edges [Hz]'};
vectLenMax = max(EEG.pacScan.resulVectLen(:));
miMax      = max(EEG.pacScan.lfoPhaseMi(:));
disp('Data are stored in EEG.pacScan.')


if plotType == 3
    disp('Plotting skipped.')
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plot mean resultant vector length / modulation index %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figCounter  = 0;
    plotCounter = 0;
    
    if EEG.nbchan < 26
        subplotRow = ceil(EEG.nbchan/5);
    else
        subplotRow = 6;
    end
    
    if EEG.nbchan < 5
        subplotColumn = EEG.nbchan;
    else
        subplotColumn = 5;
    end
    
    for chIdx = 1:EEG.nbchan
        
        if figCounter   == 0 || plotCounter == 30;
            figCounter  = figCounter+1;
            plotCounter = 0;
            myfigure(EEG, figCounter, phaseFreqRange, ampFreqRange, plotType);
        end
        
        % open subplot
        plotCounter = plotCounter+1;
        subplot(subplotRow,subplotColumn,plotCounter)
        
        % plot
        switch plotType
            case 1
                imagesc(squeeze(EEG.pacScan.resulVectLen(chIdx,:,:)))
                colormap('jet')
                switch normalizeColor
                    case 1
                        set(gca, 'CLim', [0 vectLenMax]);
                end
            case 2
                imagesc(squeeze(EEG.pacScan.lfoPhaseMi(chIdx,:,:)))
                colormap('jet')
                switch normalizeColor
                    case 1
                        set(gca, 'CLim', [0 miMax]);
                end
        end
        axis xy
        
        
        % set up x axis
        tmpLength = max(1,round((length(ampFreqList)-1)/6));
        xTick      = 1:tmpLength:length(ampFreqList)-1;
        xTickLabel = ampFreqList(xTick);
        xTickLabel = round(xTickLabel*10)/10;
        set(gca, 'XTick', xTick,  'XTickLabel', xTickLabel)
        
        % set up y axis
        tmpLength = max(1,round((length(phaseFreqList)-1)/6));
        yTick      = 1:tmpLength:length(phaseFreqList)-1;
        yTickLabel = phaseFreqList(yTick);
        yTickLabel = round(yTickLabel*10)/10;
        set(gca, 'YTick', yTick,  'YTickLabel', yTickLabel)
        
        
        if chIdx == 1
            currentPosition = get(gca, 'position');
            colorbar
            set(gca, 'position', currentPosition);
            set(get(gca, 'XLabel'), 'String', 'HFO frequency (Hz)')
            set(get(gca, 'YLabel'), 'String', 'Phase frequency (Hz)')
        end
        
        
        % title
        str = sprintf('Ch%d', chIdx);
        title(str)
    end
    
    % finalize the halfway figure
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');
end

% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% plot phase angle %%%
% %%%%%%%%%%%%%%%%%%%%%%%%
% figCounter  = 0;
% plotCounter = 0;
% 
% if EEG.nbchan < 26
%     subplotRow = ceil(EEG.nbchan/5);
% else
%     subplotRow = 6;
% end
% 
% if EEG.nbchan < 5
%     subplotColumn = EEG.nbchan;
% else
%     subplotColumn = 5;
% end
% 
% for n = 1:EEG.nbchan
%     plotType = 'Preferred phase angles';
%     if figCounter   == 0 || plotCounter == 30;
%         figCounter  = figCounter+1;
%         plotCounter = 0;
%         myfigure(EEG, figCounter, phaseFreqRange, ampFreqRange, plotType);
%     end
%         
%     % open subplot
%     plotCounter = plotCounter+1;
%     subplot(subplotRow,subplotColumn,plotCounter)
% 
%     % plot
%     imagesc(squeeze(EEG.pacScan.pacAngle (n,:,:)))
% 
%     axis xy
% 
%     % set up x axis
%     tmpLength = round(length(hasRateList)/7);
%     xTick      = 1:tmpLength:length(hasRateList);
%     xTickLabel = hasRateList(xTick);
%     xTickLabel = round(xTickLabel*10)/10;
%     set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel)
% 
%     % set up y axis
%     tmpLength = round((length(phaseFreqList)-1)/6);
%     yTick      = 1:tmpLength:length(phaseFreqList)-1;
%     yTickLabel = phaseFreqList(yTick);
%     yTickLabel = round(yTickLabel*10)/10;
%     set(gca, 'YTick', yTick,  'YTickLabel', yTickLabel)
% 
%     % title
%     str = sprintf('Ch%d', n);
%     title(str)
% end
% 
% % finalize the halfway figure
% axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');


    
function myfigure(EEG, figCounter, phaseFreqRange, ampFreqRange, plotType)
    
    % say good bye to the current figure
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');

    figure
    set(gcf, 'Name', 'pac_scanHfoByPhase()', 'NumberTitle','off')
    set(gcf,'Color', [0.93 0.96 1]);
    maxFigNum = ceil(EEG.nbchan/30);
    switch plotType
        case 1
            plotType = 'Mean vector length';
        case 2
            plotType = 'Modulation Index';
    end
    if ampFreqRange(2) == 0
        ampFreqRange(2) = round(EEG.srate/2);
    end
    annotation(gcf,'textbox', [0 1 1 0], 'String',...
        {['Phase (LFO) vs. Amplitude (HFO) ' plotType ' (Page ' num2str(figCounter) ' of ' num2str(maxFigNum) ')']}, ...
          'HorizontalAlignment','center', 'FontSize',18, 'FitBoxToText','off', 'LineStyle','none');

function out = segmentHilbert(data,cuts,valid,part)
% Analytic amplitude or phase within each boundary segment; valid columns only.
out = zeros(size(data,1),nnz(valid));
column = 0;
for s = 1:numel(cuts)-1
    span = cuts(s):cuts(s+1)-1;
    keep = valid(span);
    if ~any(keep), continue, end
    analytic = hilbert(double(data(:,span)).').';
    if strcmp(part,'abs')
        values = abs(analytic);
    else
        values = angle(analytic);
    end
    out(:,column+(1:nnz(keep))) = values(:,keep);
    column = column+nnz(keep);
end
