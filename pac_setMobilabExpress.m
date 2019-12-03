% pac_setMobilabExpress(): removes paths to mobilab and its subfolders and
%                          sets path to mobilabExpress.

% History
% 04/19/13 ver 1.1 by Makoto. propertyGrid added.
% 03/27/13 ver 1.0 by Makoto. Created.

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

function pac_setMobilabExpress

if ~isempty(which('eegplugin_mobilab'))
    eeglabPath = which('eeglab');
    path2remove = [eeglabPath(1:end-8) 'plugins/mobilab'];
    rmpathsub(path2remove);
    disp('The paths for mobilab and its subfolders are removed')
end

if isempty(which('streamBrowserNG'))
    eeglabPath = which('eeglab');
    path2add   = [eeglabPath(1:end-8) 'plugins/PACT/mobilabExpress/'];
    path2add   = [path2add ';' eeglabPath(1:end-8) 'plugins/PACT/mobilabExpress/propertyGrid/'];
    addpath(path2add);
    disp('The path for mobilabExpress under PACT is set')
end

function rmpathsub(directory)
% RMPATHSUB Remove all subdirectories from MATLAB path.
% RMPATHSUB(DIRECTORY) removes DIRECTORY and all
% its subdirectories from the MATLAB search path.

% get path as long string
p=path;

% divide string to directories, don't
% forget the first or the last...
delim=[0 strfind(p, ':') length(p)+1];

for i=2:length(delim)
    direc = p(delim(i-1)+1:delim(i)-1);
    if strncmpi(direc, directory, length(directory))
        rmpath(direc);
    end
end