% pac_vis_eventRelatedPac: Plots time-series of mean high-frequency
%                          amplitude and low-frequncy phase.
%
% Usage:
%   >>  pac_vis_eventRelatedPac(EEG);

% Author: Makoto Miyakoshi. Cincinnati Children's Hospital 
% History
% 08/14/2024 Makoto. Created. Requested by Henrico.

% Copyright (C) 2024, Makoto Miyakoshi, Henrico Stam.
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

function pac_vis_eventRelatedPac(EEG)

% Check if the precalculated data exist.
if ~isfield(EEG, 'pac')
    warndlg('Event-related PAC not calculated.')
    return
elseif ~isfield(EEG.pac, 'windowMeanAmp')
    warndlg('Event-related PAC not calculated.')
    return
end

% Launch user input dialogue window.
prompt   = {'Enter the channel index [n]'};
dlgtitle = 'Event-related PAC visualizer';
dims     = [1 30]; % The size of the input field (rows, columns)
defaultInput = {'1'}; % Default input values
chanIdxCell = inputdlg(prompt, dlgtitle, dims, defaultInput);

% Obtain the data.
chIdx = str2double(chanIdxCell{1});
meanAmp   = EEG.pac.windowMeanAmp(  chIdx,:);
meanPhase = EEG.pac.windowMeanPhase(chIdx,:);
stdAmp   = EEG.pac.windowStdAmp(  chIdx,:);
stdPhase = EEG.pac.windowStdPhase(chIdx,:);

ampEnvUpper = meanAmp+stdAmp;
ampEnvLower = meanAmp-stdAmp;
phaseEnvUpper = meanPhase+stdPhase;
phaseEnvLower = meanPhase-stdPhase;

plotTimeLength = size(meanAmp,2);
plotTime       = 1000/EEG.srate*(1:plotTimeLength);
plotTime = plotTime-plotTime(length(plotTime)/2);

figure('NumberTitle','off', 'Name', 'pac_vis_eventRelatedPac')
yyaxis left
fill([plotTime plotTime(end:-1:1)], [phaseEnvUpper phaseEnvLower(end:-1:1)], [0 0 1], 'lineStyle', 'none', 'FaceAlpha', 0.1)
hold on
phaseLineHandle = plot(plotTime, meanPhase, 'b', 'linestyle', '-', 'linewidth', 2);
xlabel('Latency (ms)')
ylabel('Phase (radians)')
set(gca, 'ytick', [-pi 0 pi], 'YTickLabel', {'-pi' '0' 'pi'})

yyaxis right
hold on
fill([plotTime plotTime(end:-1:1)], [ampEnvUpper ampEnvLower(end:-1:1)], [1 0 0], 'lineStyle', 'none', 'FaceAlpha', 0.1)
ampLineHandle = plot(plotTime, meanAmp, 'r', 'linestyle', '-', 'linewidth', 2);
ylabel('Amplitude (\muV)')

line([0 0], ylim, 'color', [0 0 0], 'linestyle', ':')

legend([phaseLineHandle ampLineHandle], {sprintf('Phase+/-1SD (%.1f-%.1fHz)', EEG.pac.lfoPhase(1), EEG.pac.lfoPhase(2)) sprintf('Amp+/-1SD (%.1f-%.1fHz)', EEG.pac.hfoAmp(1), EEG.pac.hfoAmp(2))}, 'location', 'best')

title(sprintf('Ch %d, %d trials, MI=%.2f, p=%.3f (FDR)', chIdx, size(EEG.pac.winEdgesAfterRejection,1), EEG.pac.mi(chIdx), EEG.pac.FDR.MIpval(chIdx)))