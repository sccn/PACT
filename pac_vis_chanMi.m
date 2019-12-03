% pac_vis_chanMi() - Plot modulation index that is computed from selected
%                    datapoints by specified rate of highest amplitude in
%                    high-frequency oscillation.
% Usage:
%   >> pac_vis_chanMi(EEG);

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD 2012-
% History:
% 06/20/2013 ver 2.0 by Makoto. Mean angle and Resultant Vector Length added.
% 04/16/2013 ver 1.2 by Makoto. Confidence Interval option added.
% 12/24/2012 ver 1.1 by Makoto. Minor change added.
% 10/23/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012 Makoto Miyakoshi, JSPS/SCCN,INC,UCSD
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

function pac_vis_chanMi(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualize pac across channels %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf,'Color', [0.93 0.96 1]);

colorScale = jet(360);
numChan = EEG.nbchan;
for n = 1:numChan
    tmp1 = EEG.pac.mi;
    mask = setdiff(1:numChan,n);
    tmp2 = tmp1;
    tmp2(mask) = 0;
    
    hold on
    h1 = bar(tmp2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% color code mean phase angle %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmpAngle = ceil(circ_rad2ang(EEG.pac.angleMean(n)));
    set(h1, 'facecolor', colorScale(tmpAngle,:))
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overplot confidence interval %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
switch EEG.pac.confIntType
    case 1
        bar(EEG.pac.moduInd95CI,'FaceColor','none','EdgeColor','k');
    case 2
        bar(EEG.pac.moduInd99CI,'FaceColor','none','EdgeColor','k');
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overplot signficance %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
offset = max(EEG.pac.mi)*3/100;
sigChan = find(EEG.pac.currentMIpval < EEG.pac.alpha);
plot(sigChan, EEG.pac.mi(sigChan,1)+offset, 'LineStyle','none', 'Marker', '*', 'MarkerSize', 8, 'Color', 'r')
hold off

set(get(gca, 'XLabel'), 'String', 'Channels', 'FontSize', 18)

set(gca, 'clim', [0 360], 'FontSize', 16)
h2 = colorbar('location','northoutside');
set(h2, 'XTick', [0 90 180 270 360],'FontSize', 16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% overplot resultant vector length %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
[ax,h3,h4] = plotyy(1:EEG.nbchan, zeros(1,EEG.nbchan), 1:EEG.nbchan, EEG.pac.vectLength);
set(h3, 'Color', [1 1 1], 'LineStyle', 'none')
set(ax(1), 'YTick', 5:5:100, 'YColor', 'k', 'Box', 'off')
set(get(gca, 'YLabel'), 'String', 'Modulation Index with CI', 'FontSize', 18)
xlim(ax(1), [0 EEG.nbchan+1])
ylim(ax(1), [0 max(EEG.pac.mi)*1.1])

set(h4, 'Color', [0 0 0.75], 'LineStyle', ':', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10)
set(ax(2), 'XTick',[], 'YTick', 0.2:0.2:1, 'YColor', [0 0 0.75], 'FontSize', 16)
set(get(ax(2), 'YLabel'), 'String', 'Resultant Vector Length', 'FontSize', 18)
xlim(ax(2), [0 EEG.nbchan+1])
ylim(ax(2), [0 max(EEG.pac.vectLength)*1.1])

annotation(gcf,'textbox', [0 1 1 0],...
    'String',{'Modulation Index (bar; red star, significant), Mean Angle (color-coded), and Resultant Vector Length (line).'},...
    'HorizontalAlignment','center', 'FontSize',18, 'FitBoxToText','off', 'LineStyle','none', 'FontName','Arial');

% display explanations in the figure name
set(gcf, 'Name', 'pac_vis_chanMi()', 'numbertitle','off');

% axcopy
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');