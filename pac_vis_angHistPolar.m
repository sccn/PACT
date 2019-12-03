% pac_vis_angHistPolar: Plots angular histogram in polar plot. There are two
%                       null hypotheses tested; 1. mean angle distribution is
%                       uniform (result shown as red title with *), 2. each
%                       channel has the same mean angle distribution as that
%                       of all-channel average (result shown as red background).
%
% Usage:
%   >>  pac_vis_angHistPolar(EEG);

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

function pac_vis_angHistPolar(EEG)

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
%%% plot polarplot of probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigVectLen  = find(EEG.pac.currentPhaseTest    < EEG.pac.alpha);
sigWtsnWill = find(EEG.pac.currentWtsnWillPval < EEG.pac.alpha);
for n = 1:EEG.nbchan
    
    if figCounter   == 0 || plotCounter == 30;
        figCounter  = figCounter+1;
        plotCounter = 0;
        myfigure(EEG, figCounter, EEG.pac.currentPhaseTestName, EEG.pac.currentMcompName);
    end

    plotCounter = plotCounter+1;
    subplot(subplotRow,subplotColumn,plotCounter)

    phasetmp = angle(EEG.pac.analyticEEG(n, EEG.pac.hfoIndex{n,1}));
    
    h = circ_plot(phasetmp','hist',[], EEG.pac.numPhaseBin, true,true,'linewidth',3,'color','r');


    % plot circular stat significance 1 (null hypo == flat distribution)
    if ismember(n, sigVectLen)
        str = sprintf('*Ch%d, VL%0.2f, Ang%0.0f(%0.0f), p=%0.3f', n, EEG.pac.vectLength(n), circ_rad2ang(EEG.pac.angleMean(n)), circ_rad2ang(EEG.pac.angleStd(n)), EEG.pac.currentPhaseTest(n));
        title(str, 'Color', [1 0 0])
    else
        str = sprintf('Ch%d, VL%0.2f, Ang%0.0f(%0.0f)', n, EEG.pac.vectLength(n), circ_rad2ang(EEG.pac.angleMean(n)), circ_rad2ang(EEG.pac.angleStd(n)));
        title(str)
    end
    
    % plot circular stat significance 2 (null hypo == same as mean distribution)
    if ismember(n, sigWtsnWill)
        ph = findall(h,'type','patch');
        set(ph,'facecolor',[1,0.85,0.95]);
    end
end

% finalize the halfway figure
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');



function myfigure(EEG, figCounter, phaseTestType, mCorrectType)
    
    % say good bye to the current figure
    axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');

    figure
    set(gcf, 'Name', 'pac_vis_probHist()', 'NumberTitle','off')
    set(gcf,'Color', [0.93 0.96 1]);
    maxFigNum = ceil(EEG.nbchan/30);
    annotation(gcf,'textbox', [0 1 1 0],...
        'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase. ' phaseTestType ' Test corrected with ' mCorrectType '. Red background, deveation from mean angle of all-channel mean. (page ' num2str(figCounter) ' of ' num2str(maxFigNum) ')']},...
        'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');   