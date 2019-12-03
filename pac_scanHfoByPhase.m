% pac_scanHfoByPhase() -               
% Usage:
%   >>  pac_scanHfoByPhase(EEG);

% Author: Makoto Miyakoshi JSPS/SCCN,INC,UCSD
% History:
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add folder to path %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop phase x amplitude combinations for all the channels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLoop = (length(phaseFreqList)-1)*(length(phaseFreqList)-1)*EEG.nbchan;
counter   = 0;
h = waitbar(counter/numLoop, 'Now processing...');

hfoDataLength = size(EEG.data,2)*HASnumber/100;
vectL   = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
modIdx  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);
hfoAmp  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(ampFreqList)-1);

for lfoFreqIdx = 1:length(phaseFreqList)-1
    tmpLfoPhase = filt_fir1fft(double(EEG.data'), phaseFreqList(lfoFreqIdx), phaseFreqList(lfoFreqIdx+1), EEG.srate);
    analyPhase  = angle(hilbert(tmpLfoPhase))';
    
    for hfoFreqIdx = 1:length(ampFreqList)-1
        tmpHfoAmp = filt_fir1fft(double(EEG.data'), ampFreqList(hfoFreqIdx), ampFreqList(hfoFreqIdx+1), EEG.srate);
        analyAmp  = abs(hilbert(tmpHfoAmp))';
        thresholdVector = prctile(analyAmp', 100-HASnumber);

        for chIdx = 1:EEG.nbchan
            currentChHfoMask = analyAmp(chIdx,:)>=thresholdVector(chIdx);
            tmpChHfo   = analyAmp(  chIdx,currentChHfoMask);
            tmpChPhase = analyPhase(chIdx,currentChHfoMask);
            vectL( chIdx,lfoFreqIdx,hfoFreqIdx) = circ_r(tmpChPhase');
            mAngle(chIdx,lfoFreqIdx,hfoFreqIdx) = circ_rad2ang(circ_mean(tmpChPhase'));
            modIdx(chIdx,lfoFreqIdx,hfoFreqIdx) = abs(mean(tmpChHfo.*exp(1i*tmpChPhase)));
            hfoAmp(chIdx,lfoFreqIdx,hfoFreqIdx) = mean(tmpChHfo);
            
            counter = counter+1;
            waitbar(counter/numLoop, h);
        end
    end
end
close(h)

% Store the results.
EEG.pacScan.lfoPhaseVect = single(vectL);
EEG.pacScan.pacAngle     = single(mAngle);
EEG.pacScan.lfoPhaseMi   = single(modIdx);
EEG.pacScan.meanHfoAmp   = single(hfoAmp);
EEG.pacScan.phaseFreqEdges = phaseFreqList;
EEG.pacScan.ampFreqEdges   = ampFreqList;
EEG.pacScan.dataDimensions = {'Channels', 'Phase Freq Edges [Hz]', 'HFO Freq Edges [Hz]'};
vectLenMax = max(EEG.pacScan.lfoPhaseVect(:));
miMax      = max(EEG.pacScan.lfoPhaseMi(:));
disp('Data are stored in EEG.pacScan.')

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