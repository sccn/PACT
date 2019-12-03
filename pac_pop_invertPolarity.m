% pac_pop_invertPolarity(): invert EEG polarity for clinical convention

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD
% History
% 05/01/2014 ver 1.2 by Makoto. Debug (if EEG.pac does not exist)
% 03/30/2014 ver 1.1 by Makoto. Revised minorly.
% 03/22/2013 ver 1.0 by Makoto. Created.

% Copyright (C) 2013, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
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

function EEG = pac_pop_invertPolarity(EEG)

% Create EEG.pac = [] if EEG.pac does not exist
if ~isfield(EEG, 'pac')
    EEG.pac = [];
end

% Create EEG.pac.polarityInverse == 0 if polarityInverse does not exist
if ~isfield(EEG.pac, 'polarityInverse')
    EEG.pac.polarityInverse = 1;
end

% Retrieve the current polarity state
if EEG.pac.polarityInverse==1
    STRING = 'non-inverted';
else
    STRING = 'inverted';
end

% Ask user whether or not invert polarity
userInput = questdlg(['Currently EEG is ' STRING '. Invert EEG polarity?']);

% if Yes, invert polarity
if strcmp(userInput, 'Yes')
    EEG.data = EEG.data*-1;
    EEG = pop_editset(EEG, 'setname', [EEG.setname ' inverted']);
    disp(sprintf('\nEEG polarity is inverted.'))
    EEG.pac.polarityInverse = EEG.pac.polarityInverse*-1;
else
    disp(sprintf('\nEEG polarity is NOT inverted.'))
    error('Polarity inversion canceled.')
end

return