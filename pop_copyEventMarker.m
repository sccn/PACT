% pop_copyEventMarker() - Copy event markers from the current dataset to
%                         the selected dataset. 
% Usage:
%   >> [ALLEEG,EEG] = pop_copyEventMarker(ALLEEG,EEG);

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD 2012-
% History:
% 12/26/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012 Makoto Miyakoshi, JSPS/SCCN,INC,UCSD;
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

function [ALLEEG,EEG] = pop_copyEventMarker(ALLEEG,EEG)

userInput = inputgui('title', 'pac_pop_main', 'geom', ...
   {{2 1 [0 0] [1 1]}   {2 1 [1 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Which dataset to copy current event markers?'} {'style' 'edit' 'string' ''}});

% if canceled, escape
if isempty(userInput), return, end

% decode user input
destination = str2double(userInput{1,1});

% Validate destination without mutating other datasets on invalid input.
if ~isscalar(destination) || ~isfinite(destination) || ...
        destination~=round(destination) || destination<1 || destination>numel(ALLEEG)
    error('Invalid dataset index.');
end
target = ALLEEG(destination);
if EEG.pnts ~= target.pnts || EEG.srate ~= target.srate || ...
        EEG.trials ~= target.trials || ~isequal(size(EEG.data),size(target.data))
    error('Source and destination data dimensions or sampling rates differ.');
end

% Event information and urevents must be copied as a coherent pair.
target.event = EEG.event;
target.urevent = EEG.urevent;
target = eeg_checkset(target,'eventconsistency');
ALLEEG(destination)=target;
fprintf('Event markers copied to dataset %d.\n',destination);
