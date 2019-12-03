% pac_vis_rawDataPlot(): Plots raw EEG data (indigo) overwritten with
%                        HFO-indexed data (red). 
% Usage:
%   >> EEG = pac_vis_rawDataPlot(EEG);
%
% Inputs:
%  EEG     : EEGLAB structure

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD 2012-
% History
% 04/16/2013 ver 1.3 by Makoto. high-pass filter bug fixed.
% 12/24/2012 ver 1.2 by Makoto. Minor change added.
% 10/22/2012 ver 1.1 by Makoto.
% 09/24/2012 ver 1.0 by Makoto.

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

function pac_vis_rawDataPlot(EEG)

userInput = inputgui('title', 'pac_vis_rawDataPlot()', 'geom', ...
   {{2 1 [0 0] [1 1]}   {2 1 [1 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Apply HFO filter?'} {'style' 'popupmenu' 'string' 'No|Yes' 'tag' 'filter' 'value' 1}});

if userInput{1,1} == 2
    EEG.data = abs(EEG.pac.analyticEEG);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create an initial event %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(EEG.event)
    EEG.event.type       = 'start';
    EEG.event.latency    = 1;
    EEG.event.init_index = 1;
    EEG.event.init_time  = 1/EEG.srate;

    EEG.urevent.type       = 'start';
    EEG.urevent.latency    = 1;
    EEG.urevent.init_index = 1;
    EEG.urevent.init_time  = 1/EEG.srate;
end

%%%%%%%%%%%%%%%%%%%%%%
%%% create redMask %%%
%%%%%%%%%%%%%%%%%%%%%%
redMask = EEG.data;
for n = 1:EEG.nbchan
    nanMask = setdiff(1:length(EEG.data), EEG.pac.hfoIndex{n,1});
    redMask(n,nanMask) = NaN;
end

%%%%%%%%%%%%
%%% plot %%%
%%%%%%%%%%%%
% pop_eegplot(EEG.data, [], 0, 0, [], 'data2', redMask);
eegplot(EEG.data, 'srate', EEG.srate, 'data2', redMask)