% pac_setPath() - set up path for PACT.

% Author:
% Makoto Miyakoshi, 2014 Swartz Center for Computational Neuroscience, INC, UCSD
%
% History:
% 03/28/2014 ver 1.00 by Makoto. Created.

% Copyright (C) 2014, Makoto Miyakoshi SCCN, INC, UCSD
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function ok = pac_setPath

% set the check flag false
ok = false;

% find the path to itself
pactRoot = fileparts(which('pac_setPath.m'));
        
% add subfolders to path (if not already present)
if ~exist('32px-Dialog-apply.svg.png','file')
    addpath(genpath(pactRoot));
end

ok = true;