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

if isfield(EEG, 'pacScan')
    EEG = rmfield(EEG,'pacScan');
end

% Prepare log-scaled frequency edges.
phaseFreqList = logspace(log10(phaseFreqRange(1)), log10(phaseFreqRange(2)), numFreqBins);
ampFreqList   = logspace(log10(ampFreqRange(1)),   log10(ampFreqRange(2)),   numFreqBins);
phaseFreqTbw  = mean(diff(phaseFreqList))/2; % /2 is a heuristic solution to improve freq resolution. May need to be removed if it makes the filter too aggressive.
ampFreqTbw    = mean(diff(ampFreqList));

% Determine filter parametes for phase. (taken from pop_firwsord.m)
df  = phaseFreqTbw/EEG.srate; % Normalize transition band width
wtype = 4; % wtypes = {'rectangular' 'bartlett' 'hann' 'hamming' 'blackman' 'kaiser'};
devs = {0.089 0.056 0.0063 0.0022 0.0002 []};
dev = devs{wtype};
dfs = [0.9 2.9 3.1 3.3 5.5]; 
m = dfs(wtype) / df;
phaseFilterOrder = ceil(m / 2) * 2;
disp(sprintf('For phase analysis, FIR window type, Hamming, TBW, %.3f.', phaseFreqTbw))
disp(sprintf('Band-pass filter center frequencies are shown below.'))
phaseFreqList

% Determine filter parametes for amp.
df  = ampFreqTbw/EEG.srate; % Normalize transition band width
wtype = 4; % wtypes = {'rectangular' 'bartlett' 'hann' 'hamming' 'blackman' 'kaiser'};
devs = {0.089 0.056 0.0063 0.0022 0.0002 []};
dev = devs{wtype};
dfs = [0.9 2.9 3.1 3.3 5.5]; 
m = dfs(wtype) / df;
ampFilterOrder = ceil(m / 2) * 2;
disp(sprintf('For Amp analysis, FIR window type, Hamming, TBW, %.3f.', ampFreqTbw))
disp(sprintf('Band-pass filter center frequencies are shown below.'))
ampFreqList



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove 1-s peri-boundary data points. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is mainly for window-rejected data to avoid filter's edge effect (11/04/2020 Makoto. Added).
if ~isempty(EEG.event)
    boundaryIdx = find(strcmp({EEG.event.type}, 'boundary'));
    if ~isempty(boundaryIdx)
        boundaryLatency = round([EEG.event(boundaryIdx).latency]);
        periBoundaryIdx = [1:EEG.srate (1+EEG.pnts-EEG.srate):EEG.pnts]; % The first and the last 1-s by the edge effect.
        for boundaryIdx = 1:length(boundaryLatency)
            periBoundaryIdx = [periBoundaryIdx 1+boundaryLatency(boundaryIdx)-EEG.srate:boundaryLatency(boundaryIdx)+EEG.srate];
        end
        periBoundaryIdx = unique(periBoundaryIdx);
        
        % Trim indices outside the data range.
        periBoundaryIdx(periBoundaryIdx<1) = [];
        periBoundaryIdx(periBoundaryIdx>EEG.pnts) = [];
    else
        boundaryLatency = [];
        periBoundaryIdx = [1:EEG.srate (1+EEG.pnts-EEG.srate):EEG.pnts]; % The first and the last 1-s by the edge effect.
    end
else
    boundaryLatency = [];
    periBoundaryIdx = [1:EEG.srate (1+EEG.pnts-EEG.srate):EEG.pnts]; % The first and the last 1-s by the edge effect.
end
disp(sprintf('%.0f boundaries found. %.1f s data will be removed.', length(boundaryLatency), length(periBoundaryIdx)/EEG.srate))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop phase x amplitude combinations for all the channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLoop       = (length(phaseFreqList)-1)*(length(phaseFreqList)-1)*EEG.nbchan;
counter       = 0;
% waitbarHandle = waitbar(counter/numLoop, 'Now processing...');

hfoDataLength = size(EEG.data,2)*HASnumber/100;
resulVectLen   = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
mAngle  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
modIdx  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
hfoAmp  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
    
EEGorig = EEG;
for phaseFreqIdx = 1:length(phaseFreqList)-1
    
    % Apply band-pass filter to extract low-freq phase.
    EEG = pop_firws(EEGorig, 'fcutoff', [phaseFreqList(phaseFreqIdx), phaseFreqList(phaseFreqIdx+1)], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', phaseFilterOrder, 'minphase', 0);
    analyPhase = angle(hilbert(EEG.data'))'; % hilbert() works columnwise!!! (11/23/2020)
    analyPhase(:,periBoundaryIdx) = [];
    
    % Apply band-pass filter to extract high-freq amp.
    for ampFreqIdx = 1:length(ampFreqList)-1
        
        EEG = pop_firws(EEGorig, 'fcutoff', [ampFreqList(ampFreqIdx), ampFreqList(ampFreqIdx+1)], 'ftype', 'bandpass', 'wtype', 'hamming', 'forder', ampFilterOrder, 'minphase', 0);
        analyAmp  = abs(hilbert(EEG.data'))'; % hilbert() works columnwise!!! (11/23/2020)
        analyAmp(:,periBoundaryIdx) = [];
        
        thresholdVector = prctile(analyAmp', 100-HASnumber);

        for chIdx = 1:EEG.nbchan
            currentChHfoMask = analyAmp(chIdx,:)>thresholdVector(chIdx);
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
            
            [pval, zStats] = circ_rtest(tmpChPhase');
            rayleighZ(   chIdx,phaseFreqIdx,ampFreqIdx) = zStats;
            resulVectLen(chIdx,phaseFreqIdx,ampFreqIdx) = circ_r(tmpChPhase');
            mAngle(      chIdx,phaseFreqIdx,ampFreqIdx) = circ_rad2ang(circ_mean(tmpChPhase'));
            modIdx(      chIdx,phaseFreqIdx,ampFreqIdx) = abs(mean(tmpChHfo.*exp(1i*tmpChPhase)));
            hfoAmp(      chIdx,phaseFreqIdx,ampFreqIdx) = mean(tmpChHfo);
            
            % To store analytic phase for Daisuke's paper. (11/20/2020)
            if ampFreqIdx == 1 & chIdx == 1
                    analyticPhaseStack(chIdx, phaseFreqIdx, ampFreqIdx,:) = tmpChPhase;
            end
            
            if length(tmpChPhase) == size(analyticPhaseStack,4)
                analyticPhaseStack(chIdx, phaseFreqIdx, ampFreqIdx,:) = tmpChPhase;
            elseif length(tmpChPhase) > size(analyticPhaseStack,4)
                analyticPhaseStack(chIdx, phaseFreqIdx, ampFreqIdx, :) = tmpChPhase(1:size(analyticPhaseStack,4));
            elseif length(tmpChPhase) < size(analyticPhaseStack,4)
                lengthDifference = size(analyticPhaseStack,4)-length(tmpChPhase);
                analyticPhaseStack(chIdx, phaseFreqIdx, ampFreqIdx, :) = [tmpChPhase tmpChPhase(1:lengthDifference)];
            end

            counter = counter+1;
%             waitbar(counter/numLoop, waitbarHandle);
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
EEG.pacScan.analyticPhaseStack = single(analyticPhaseStack);
EEG.pacScan.phaseFreqEdges = phaseFreqList;
EEG.pacScan.ampFreqEdges   = ampFreqList;
EEG.pacScan.phaseFreqTbw   = phaseFreqTbw;
EEG.pacScan.ampFreqTbw     = ampFreqTbw;
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
                imagesc(squeeze(EEG.pacScan.lfoPhaseVect(chIdx,:,:)))
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
        tmpLength = round((length(ampFreqList)-1)/6);
        xTick      = 1:tmpLength:length(ampFreqList)-1;
        xTickLabel = ampFreqList(xTick);
        xTickLabel = round(xTickLabel*10)/10;
        set(gca, 'XTick', xTick,  'XTickLabel', xTickLabel)
        
        % set up y axis
        tmpLength = round((length(phaseFreqList)-1)/6);
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