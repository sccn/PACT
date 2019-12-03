% pac_pop_countHfo()

% Author: Makoto Miyakoshi. SCCN,INC,UCSD
% History:
% 03/31/2014 ver 1.0 by Makoto. Created.

% Copyright (C) 2014 Makoto Miyakoshi, SCCN,INC,UCSD
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

function pac_pop_countHfo(EEG)

% if EEG.event does not any event
if isempty(EEG.event)
    error('No event.')
    return
end

% extract all event types
allEventTypes = {EEG.event.type}';

% extract all event types that says 'HFO'
hfoMatchIdx = find(strcmp(allEventTypes,'HFO'));

% if no HFO
if ~any(hfoMatchIdx)
    error('No type ''HFO''.')
    return
end

% extract HFO marking frequency in all channels
allHfoChann = cell2mat({EEG.event(hfoMatchIdx).channel}');

% count how many HFO markers each channel has
hfoCountList = zeros(EEG.nbchan,1);
for ch = 1:EEG.nbchan
    hfoCountList(ch) = sum(allHfoChann==ch);
    disp(sprintf('Channel %3.0f, %3.0f HFOs marked', ch, hfoCountList(ch)))
end

% export as .txt
[filename,filepath] = uiputfile('*.txt', 'Enter save file name -- pac_pop_countHfo()'); 
drawnow
if filename == 0
    disp('Ascii file export canceled.')
    return;
end
filename = [filepath filename];
save(filename, 'hfoCountList', '-ascii');
disp(['Exported to the following folder: ' filename])