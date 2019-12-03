% pac_vis_phaseSortAmp: Plots phase-bin-mean amplitudes. 
%
% Usage:
%   >>  pac_vis_phaseSortAmp(EEG);

% Author: Makoto Miyakoshi, Arnaud Delorme JSPS/SCCN,INC,UCSD
% History
% 01/31/2013 ver 1.2 by Makoto. yaxis scale adjusted.
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

function pac_vis_phaseSortAmp(EEG)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot phase-sorted amplitude %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for n = 1:EEG.nbchan
    if figCounter   == 0 | plotCounter == 30;
        figCounter  = figCounter+1;
        plotCounter = 0;
        myfigure(EEG, figCounter, EEG.pac.currentAmpTestName, EEG.pac.currentMcompName)
    end
    plotCounter = plotCounter+1;
    subplot(subplotRow,subplotColumn,plotCounter)

    bar(EEG.pac.phaseSortedAmp{n,1}, 'FaceColor', [0 86/255 69/255], 'BarWidth', 1); % Billiard 
    hold on
    errorbar(1:EEG.pac.numPhaseBin, EEG.pac.phaseSortedAmp{n,1}, zeros(1,EEG.pac.numPhaseBin), EEG.pac.phaseSortedAmpSe{n,1}, 'LineStyle','none','Color',[0 0 0])
    
    yMax = max(EEG.pac.phaseSortedAmp{n,1} + EEG.pac.phaseSortedAmpSe{n,1});
    yMin = min(EEG.pac.phaseSortedAmp{n,1});
    xlim([0.5 EEG.pac.numPhaseBin+0.5])
    ylim([yMin*0.85 yMax*1.15])
    xlabel('phase')
    ylabel('amp (uV)')
    hold on

    % overplot a cosine wave
    t = 0:0.1:EEG.pac.numPhaseBin;
    oneCycleConstant = EEG.pac.numPhaseBin/6.28; % 0 to 2pi
    plot(t, (cos(t/oneCycleConstant)+1)/4*(yMax*1.15-yMin*0.85) + yMin*0.85, 'r', 'LineWidth',2)
    set(gca, 'XTickLabel',{'0','pi','2*pi'},'XTick',[0 EEG.pac.numPhaseBin/2 EEG.pac.numPhaseBin], 'CLim',[1 2])

    % plot statistical significance
    sigChan = find(EEG.pac.currentAmpTest < EEG.pac.alpha);
    if ismember(n, sigChan)
        str = sprintf('*Ch%d, p=%0.3f', n, EEG.pac.currentAmpTest(n));
        title(str, 'Color', [1 0 0])
    else
        str = sprintf('Ch%d', n);
        title(str)
    end
end

% say good bye to the current figure
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');


function myfigure(EEG, figCounter, distribTestType, mCorrectType)
    
    % say good bye to the current figure
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');

    figure
    set(gcf, 'Name', 'pac_vis_phaseSortAmp()', 'NumberTitle','off')
    set(gcf,'Color', [0.93 0.96 1]);
    maxFigNum = ceil(EEG.nbchan/30);
    annotation(gcf,'textbox', [0 1 1 0],...
        'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase. ' distribTestType ' Test corrected with ' mCorrectType '.  (page ' num2str(figCounter) ' of ' num2str(maxFigNum) ')']},...
        'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');  