% pac_pop_compareLfoFreqs: The wrapper for pac_compareLfoFreqs(). Launches
%                          GUI to collect user input and pass them to the
%                          main funciton.
%
% Usage:
%   >>  pac_pop_compareLfoFreqs(EEG);

% Author: Makoto Miyakoshi SCCN,INC,UCSD
% History
% 06/24/2013 ver 1.0 by Makoto. Created.

% Copyright (C) 2013, Makoto Miyakoshi SCCN,INC,UCSD
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

function pac_pop_compareLfoFreqs(EEG)

userInput = inputgui('title', 'pac_pop_compareLfoFreqs', 'geom', ...
   {{2 4 [0 0] [1 1]} {2 4 [1 0] [1 1]} ...
    {2 4 [0 1] [1 1]} {2 4 [1 1] [1 1]} ...
    {2 4 [0 2] [1 1]} {2 4 [1 2] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Phase (LFO) freq range [loHz hiHz]'}             {'style' 'text' 'string' '[0.5 1 2 4 8]'}... 
    {'style' 'text' 'string' 'HFO filter [loHz hiHz]; [loHz 0] -> high-pass' } {'style' 'edit' 'string' '100 0'}...
    {'style' 'text' 'string' 'HFO amp percentile [%]'}                         {'style' 'edit' 'string' '2'}});

lfoFreqRange     = [0.5 1 2 4 8];
hfoFilter        = str2double(userInput{1,1});
hfoAmpPercentile = str2double(userInput{1,2});
plotFlag         = 1;

pac_compareLfoFreqs(EEG, lfoFreqRange, hfoFilter, hfoAmpPercentile, plotFlag);
