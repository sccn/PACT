% pac_vis_angHistBar: Plots angular histogram in bar plot. There are two
%                     null hypotheses tested; 1. mean angle distribution is
%                     uniform (result shown as red title with *), 2. each
%                     channel has the same mean angle distribution as that
%                     of all-channel average (result shown as red background).
%
% Usage:
%   >>  pac_vis_angHistBar(EEG);

% Author: Makoto Miyakoshi, Arnaud Delorme JSPS/SCCN,INC,UCSD
% History
% 06/21/2013 ver 2.0 by Makoto. SD added. Watson-Williams test added.
% 12/24/2012 ver 1.1 by Makoto. Minor change added.
% 11/01/2012 ver 1.0 by Makoto. Created.

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

function pac_vis_angHistBar(EEG)

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot histogram of probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigVectLen  = find(EEG.pac.currentPhaseTest    < EEG.pac.alpha);
sigWtsnWill = find(EEG.pac.currentWtsnWillPval < EEG.pac.alpha);
for n = 1:EEG.nbchan
    
    if figCounter   == 0 || plotCounter == 30;
        figCounter  = figCounter+1;
        plotCounter = 0;
        myfigure(EEG, figCounter, EEG.pac.currentPhaseTestName, EEG.pac.currentMcompName);
    end
    
    % open subplot
    plotCounter = plotCounter+1;
    subplot(subplotRow,subplotColumn,plotCounter)
        
    % angular histogram
    phasetmp     = angle(EEG.pac.analyticEEG(n, EEG.pac.hfoIndex{n,1}));
    phasetmpPosi = phasetmp(phasetmp>0);
    phasetmpNega = phasetmp(phasetmp<0);
    phasetmpNega = 2*pi+phasetmpNega;
    phasetmp     = [phasetmpPosi phasetmpNega];
    hist(phasetmp, EEG.pac.numPhaseBin);
    h = findobj(gca, 'Type', 'patch');
    set(h, 'FaceColor', [16/255 67/255 100/255]); % Marine Blue
    hold on;

    % overplot a cosine wave
    t = 0:0.01:2*pi;
    plot(t, (cos(t)+1)/4*max(EEG.pac.histBinMax), 'r', 'LineWidth',2)
    set(gca, 'XTickLabel',{'0','pi','2*pi'},'XTick',[0 3.1415 6.2832], 'CLim',[1 2])
    xlim([0 2*pi+0.01])
    ylim([0 max(EEG.pac.histBinMax)*1.01])
    xlabel('phase')
    ylabel(['probability*' num2str(length(phasetmp))])
    
    % plot circular stat significance 1 (null hypo == flat distribution)
    if ismember(n, sigVectLen)
        str = sprintf('*Ch%d, VL%0.2f, Rad%0.2f(%0.2f), p=%0.3f', n, EEG.pac.vectLength(n), EEG.pac.angleMean(n), EEG.pac.angleStd(n), EEG.pac.currentPhaseTest(n));
        title(str, 'Color', [1 0 0])
    else
        str = sprintf('Ch%d, VL%0.2f, Rad%0.2f(%0.2f)', n, EEG.pac.vectLength(n), EEG.pac.angleMean(n), EEG.pac.angleStd(n));
        title(str)
    end
    
    % plot circular stat significance 2 (null hypo == same as mean distribution)
    if ismember(n, sigWtsnWill)
        set(gca, 'Color', [1 0.85 0.95]);
    end
end

% finalize the halfway figure
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');


    
    
function myfigure(EEG, figCounter, phaseTestType, mCorrectType)
    
    % say good bye to the current figure
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');

    figure
    set(gcf, 'Name', 'pac_vis_angHistBar()', 'NumberTitle','off')
    set(gcf,'Color', [0.93 0.96 1]);
    maxFigNum = ceil(EEG.nbchan/30);
    annotation(gcf,'textbox', [0 1 1 0],...
        'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase. ' phaseTestType ' Test corrected with ' mCorrectType '. Red background, deveation from mean angle of all-channel mean. (page ' num2str(figCounter) ' of ' num2str(maxFigNum) ')']},...
        'HorizontalAlignment','center', 'FontSize',16, 'FitBoxToText','off', 'LineStyle','none');   
    