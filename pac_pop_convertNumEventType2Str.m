% pac_pop_convertNumelEventType2String() - Convert numeric event names (i.e. 'type') to string

% Author: Makoto Miyakoshi, SCCN,INC,UCSD
% History:
% 08/11/2014 ver 1.0 by Makoto.. 

% Copyright (C) 2014 Makoto Miyakoshi, SCCN,INC,UCSD;
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

function EEG = pac_pop_convertNumelEventType2String(EEG)

for n = 1:length(EEG.event)
    if isnumeric(EEG.event(n).type)
        userInput = questdlg('Numeric event type is detected: convert to string?');
        if     strcmp(userInput, 'Yes')
            for m = 1:length(EEG.event)
                EEG.event(m).type = num2str(EEG.event(m).type);
                disp([EEG.event(m).type ' converted.'])
            end
            disp('Done.')
            return
        else
            disp('Conversion cancelled.')
            return
        end
    end
end
disp('All event types are string.')