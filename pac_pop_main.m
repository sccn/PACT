% pac_pop_main(): Launch GUI to collect user inputs for computing
%                 phase-amplitude coupling. Type 'help pac_man' for
%                 details in computation.
% Usage:
%   >>  [ALLEEG EEG] = pac_pop_main(ALLEEG, EEG, CURRENTSET);

% Author: Makoto Miyakoshi, Arnaud Delorme JSPS/SCCN,INC,UCSD
%
% History:
% 08/09/2024 Makoto. Fix request by Henrico.
% 01/13/2021 Makoto. Checked for moving to Github.
% 09/20/2019 Makoto. eegh() supported. Thank you Brian Kavanaugh for using PACT.
% 06/19/2014 ver 2.2 by Makoto. dropDownStrings bug fixed.
% 01/09/2013 ver 2.1 by Makoto. isfield(EEG.event, 'type') added.
% 12/25/2012 ver 2.0 by Makoto. MobiLab and VisEd supported. assignin used.
% 10/19/2012 ver 1.0 by Makoto. Created.

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

function [ALLEEG, EEG] = pac_pop_main(ALLEEG, EEG, CURRENTSET)

if nargin < 3
    help pac_pop_main;
    return;
end;

% add folder to path
% -----------------------
if ~exist('circ_rtest','file')
    p = which('eegplugin_pac');
    p = p(1:strfind(p,'eegplugin_pac.m')-1);
    p = genpath(p);
    addpath(p);
end

% if EEG.event has chan (character) generated by VisEd, convert them to channel (double)
if isfield(EEG.event, 'chan')
    for n = 1:length(EEG.event)
        if ~isempty(EEG.event(1,n).chan)
            hitIndex = strcmp({EEG.chanlocs.labels}, EEG.event(1,n).chan);
            if any(hitIndex)
                EEG.event(1,n).channel = find(hitIndex);
            end
        end
    end
end

% create a letter string for drop down menu for selecting event 
if isfield(EEG.event, 'type')
    eventTypes = unique({EEG.event.type});
    for n = 1:length(eventTypes)
        if n == 1
            dropDownStrings = eventTypes{1,n};
        else
            dropDownStrings = [dropDownStrings '|' eventTypes{1,n}]; %#ok<*AGROW>
        end
    end
else % if EEG.event does not exist
    dropDownStrings = '(non-existent)';
end

try % try to find stored parameters
    userInput = inputgui('title', 'pac_pop_main()', 'geom', ...
       {{2 8 [0 0] [1 1]}   {2 8 [1 0] [1 1]} ...
        {2 8 [0 1] [1 1]}   {2 8 [1 1] [1 1]} ...
        {2 8 [0 2] [1 1]}   {2 8 [1 2] [1 1]} ...
        {2 8 [0 3] [1 1]}   {2 8 [1 3] [1 1]} ...
        {2 8 [0 4] [1 1]}   {4 8 [2 4] [1 1]} {4 8 [3 4] [1 1]} ...
        {2 8 [0 5] [1 1]}   {2 8 [1 5] [1 1]} ...
        {2 8 [0 6] [1 1]}   {2 8 [1 6] [1 1]} ...
        {2 8 [0 7] [1 1]}   {2 8 [1 7] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' 'Phase freq range [lohz hihz]'}                    {'style' 'edit' 'string' num2str(EEG.pac.lfoPhase)   } ...
        {'style' 'text' 'string' 'Amp freq range [lohz hihz]' }                     {'style' 'edit' 'string' num2str(EEG.pac.hfoAmp)     } ...
        {'style' 'text' 'string' 'HFO amplitude percentile [%]'}                    {'style' 'edit' 'string' num2str(EEG.pac.hfoTopRatio)} ...
        {'style' 'text' 'string' 'Sampling pool'}                                   {'style' 'popupmenu' 'string' 'each channel|all channels|handpick (specialized single-channel events)|handpick (EEGLAB''s events)' 'tag' 'hfoPool' 'value' EEG.pac.hfoPool} ...
        {'style' 'text' 'string' 'If handpicked, event type and win size [+/- ms]'} {'style' 'popupmenu' 'string' dropDownStrings 'value' EEG.pac.whichMarker} {'style' 'edit' 'string' num2str(EEG.pac.windowLength)} ...
        {'style' 'text' 'string' 'Significance threshold [p]'}                      {'style' 'edit' 'string' num2str(EEG.pac.alpha)      } ...
        {'style' 'text' 'string' 'Number of surrogation [N]'}                       {'style' 'edit' 'string' num2str(EEG.pac.numSurro)   } ...
        {'style' 'text' 'string' 'Number of phase bins [N]'}                        {'style' 'edit' 'string' num2str(EEG.pac.numPhaseBin)}});

catch % if the stored parameters were not found, use default 
    userInput = inputgui('title', 'pac_pop_main', 'geom', ...
       {{2 8 [0 0] [1 1]}   {2 8 [1 0] [1 1]} ...
        {2 8 [0 1] [1 1]}   {2 8 [1 1] [1 1]} ...
        {2 8 [0 2] [1 1]}   {2 8 [1 2] [1 1]} ...
        {2 8 [0 3] [1 1]}   {2 8 [1 3] [1 1]} ...
        {2 8 [0 4] [1 1]}   {4 8 [2 4] [1 1]} {4 8 [3 4] [1 1]} ...
        {2 8 [0 5] [1 1]}   {2 8 [1 5] [1 1]} ...
        {2 8 [0 6] [1 1]}   {2 8 [1 6] [1 1]} ...
        {2 8 [0 7] [1 1]}   {2 8 [1 7] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' 'Phase freqs [lohz hihz]'}                         {'style' 'edit' 'string' '1 2'   } ...
        {'style' 'text' 'string' 'Amplitude freqs [lohz hihz]; [lohz 0]->hipass' }  {'style' 'edit' 'string' '100 0'} ...
        {'style' 'text' 'string' 'HFO amplitude percentile [%]'}                    {'style' 'edit' 'string' '2'    } ...
        {'style' 'text' 'string' 'Sampling pool'}                                   {'style' 'popupmenu' 'string' 'each channel|all channels|handpicked' 'tag' 'hfoPool' 'value' 1} ...
        {'style' 'text' 'string' 'If handpicked, event type and win size [+/- ms]'} {'style' 'popupmenu' 'string' dropDownStrings 'value' 1} {'style' 'edit' 'string' ''} ...
        {'style' 'text' 'string' 'Significance threshold [p]'}                      {'style' 'edit' 'string' '0.01'  } ...
        {'style' 'text' 'string' 'Number of surrogation [N]'}                       {'style' 'edit' 'string' '2000'  } ...
        {'style' 'text' 'string' 'Number of phase bins [N]'}                        {'style' 'edit' 'string' '36'    }});
end

% if canceled, escape
if isempty(userInput), return, end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% store these parameters to EEG.pac %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG.pac.lfoPhase    = str2num(userInput{1}); %#ok<*ST2NM>
EEG.pac.hfoAmp      = str2num(userInput{2});
EEG.pac.hfoTopRatio = str2num(userInput{3});
EEG.pac.hfoPool     = userInput{4};
EEG.pac.whichMarker = userInput{5};
EEG.pac.windowLength= str2num(userInput{6});
EEG.pac.alpha       = str2num(userInput{7});
EEG.pac.numSurro    = str2num(userInput{8});
EEG.pac.numPhaseBin = str2num(userInput{9});

assignin('base','EEG',EEG);

% prepare varargin for std_envtopo
options = '';
    options = [ options '''lfoPhase'',    [' userInput{1} '],' ];
    options = [ options '''hfoAmp'',      [' userInput{2} '],' ];
    options = [ options '''hfoTopRatio'', [' userInput{3} '],' ];
    options = [ options '''hfoPool'',     [' num2str(userInput{4}) '],' ];
    options = [ options '''whichMarker'', [' num2str(userInput{5}) '],' ];
    options = [ options '''windowLength'',[' userInput{6} '],' ];
    options = [ options '''alpha'',       [' userInput{7} '],' ];
    options = [ options '''numSurro'',    [' userInput{8} '],' ];
    options = [ options '''numPhaseBin'', [' userInput{9} '],' ];
options = eval([ '{' options '}' ]);    
    
% run
EEG = pac_man(EEG, options{:});
EEG = pac_pop_statsSetUp(EEG);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% Update eegh. (09/20/2019 updated).
optionStrings = sprintf('''lfoPhase'', [%s], ''hfoAmp'', [%s], ''hfoTopRatio'', %d, ''hfoPool'', %d, ''whichMarker'', %d, ''windowLength'', [%d], ''alpha'', [%d], ''numSurro'', %d, ''numPhaseBin'', %d', ...
                           num2str(options{2}), num2str(options{4}), options{6}, options{8}, options{10}, options{12}, options{14}, options{16}, options{18});
com = sprintf('EEG = pac_man(EEG, %s)', optionStrings);
EEG = eegh(com, EEG);