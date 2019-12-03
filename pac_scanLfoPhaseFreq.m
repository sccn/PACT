% pac_scanLfoPhaseFreq() -               
% Usage:
%   >>  pac_scanLfoPhaseFreq(EEG);

% Author: Makoto Miyakoshi JSPS/SCCN,INC,UCSD
% History:
% 06/21/2013 ver 3.0 by Makoto. Filter renewed thanks to Christian's advice. Now works 10-20 times faster.
% 01/15/2013 ver 2.0 by Makoto. Renewed.
% 11/04/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
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

function EEG = pac_scanLfoPhaseFreq(EEG, phaseFreqRange, numPhaseFreqs, ampFreqRange, hasRateRange, numHasRates, plotType, normalizeColor)

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

phaseFreqList = logspace(log10(phaseFreqRange(1)), log10(phaseFreqRange(2)), numPhaseFreqs);
hasRateList   = logspace(log10(hasRateRange(1)),   log10(hasRateRange(2)),   numHasRates);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute amplitudes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nnz(ampFreqRange(2)) % band-pass filter for HFO
    filtAmp = filt_fir1fft(double(EEG.data'), ampFreqRange(1), ampFreqRange(2), EEG.srate);
else                    % high-pass filter for HFO
    filtAmp = filt_fir1fft(double(EEG.data'), ampFreqRange(1), 0, EEG.srate);
end
analyAmp             = abs(hilbert(filtAmp))';
[analyAmpSort,index] = sort(analyAmp, 2, 'descend');
edges                = ceil(hasRateList/100*length(analyAmp));
edges                = edges(end:-1:1);
initEdgeIdx          = index(:,1:edges(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop for LFO freqs, HFO percentiles, and Channels  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLoop = (length(phaseFreqList)-1)*length(edges)*EEG.nbchan;
counter   = 0;
h = waitbar(counter/numLoop, 'Now processing...');
vectL   = zeros(EEG.nbchan,length(phaseFreqList)-1,length(edges));
modIdx  = zeros(EEG.nbchan,length(phaseFreqList)-1,length(edges));

for lfoIdx = 1:length(phaseFreqList)-1
    tmpLfoPhase = filt_fir1fft(double(EEG.data'), phaseFreqList(lfoIdx), phaseFreqList(lfoIdx+1), EEG.srate);
    analyPhase  = angle(hilbert(tmpLfoPhase))';
    
    for numEdge    = 1:length(edges)
        tmpEdge    = edges(numEdge);
        tmpEdgeIdx = initEdgeIdx(:,1:tmpEdge);

        for ch = 1:EEG.nbchan
            tmpChHfo   = analyAmp(  ch,tmpEdgeIdx(ch,:));
            tmpChPhase = analyPhase(ch,tmpEdgeIdx(ch,:));
            vectL( ch,lfoIdx,length(edges)-numEdge+1) = circ_r(tmpChPhase');
            mAngle(ch,lfoIdx,length(edges)-numEdge+1) = circ_rad2ang(circ_mean(tmpChPhase'));
            modIdx(ch,lfoIdx,length(edges)-numEdge+1) = abs(mean(tmpChHfo.*exp(1i*tmpChPhase)));
            counter = counter+1;
            waitbar(counter/numLoop, h);
        end
    end
end
close(h)
            
EEG.pacScan.lfoPhaseVect = vectL ;
EEG.pacScan.pacAngle     = mAngle;
EEG.pacScan.lfoPhaseMi   = modIdx;
vectLenMax = max(EEG.pacScan.lfoPhaseVect(:));
miMax      = max(EEG.pacScan.lfoPhaseMi(:));

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

for n = 1:EEG.nbchan
    
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
            imagesc(squeeze(EEG.pacScan.lfoPhaseVect(n,:,:)))
            switch normalizeColor
                case 1
                    set(gca, 'CLim', [0 vectLenMax]);
            end
        case 2
            imagesc(squeeze(EEG.pacScan.lfoPhaseMi(n,:,:)))
            switch normalizeColor
                case 1            
                    set(gca, 'CLim', [0 miMax]);
            end
    end
    axis xy
    

    % set up x axis
    tmpLength = round(length(hasRateList)/7);
    xTick      = 1:tmpLength:length(hasRateList);
    xTickLabel = hasRateList(xTick);
    xTickLabel = round(xTickLabel*10)/10;
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel)

    % set up y axis
    tmpLength = round((length(phaseFreqList)-1)/6);
    yTick      = 1:tmpLength:length(phaseFreqList)-1;
    yTickLabel = phaseFreqList(yTick);
    yTickLabel = round(yTickLabel*10)/10;
    set(gca, 'YTick', yTick,  'YTickLabel', yTickLabel)
    
    if n == 1
        set(get(gca, 'XLabel'), 'String', 'HFO percentile (%)')
        set(get(gca, 'YLabel'), 'String', 'LFO frequency (Hz)')
    end

        
    % title
    str = sprintf('Ch%d', n);
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
    set(gcf, 'Name', 'pac_scanLfoPhaseFreq()', 'NumberTitle','off')
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
        {['HFO (' num2str(ampFreqRange(1)) '-' num2str(ampFreqRange(2)) ' Hz) coupled with ranges of LFOs. ' plotType ' is plotted. (page ' num2str(figCounter) ' of ' num2str(maxFigNum) ')']}, ...
          'HorizontalAlignment','center', 'FontSize',18, 'FitBoxToText','off', 'LineStyle','none');