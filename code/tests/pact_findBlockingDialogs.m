function offenders = pact_findBlockingDialogs(callback, root)
%PACT_FINDBLOCKINGDIALOGS PACT functions reachable from a callback that block.
%   offenders = pact_findBlockingDialogs(callbackString, rootFolder)
% Starting from the pac_*/pop_* functions named in callbackString, follows
% calls to further pac_*/pop_* functions whose files lie under rootFolder
% (functions in rootFolder are resolved there first; others via WHICH and
% ignored unless under rootFolder). Whole-line comments are ignored. Returns the names
% of functions that call inputdlg, questdlg, listdlg, uiwait or keyboard.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

blocking='\<(inputdlg|questdlg|listdlg|uiwait|keyboard)\s*\(';
callPattern='\<(pac_\w+|pop_\w+)(?=\s*\()';
pending=regexp(callback,callPattern,'match');
visited={};
offenders={};
while ~isempty(pending)
    name=pending{1};
    pending(1)=[];
    if any(strcmp(visited,name)), continue, end
    visited{end+1}=name; %#ok<AGROW>
    file=fullfile(root,[name '.m']);
    if ~isfile(file)
        file=which(name);
        if isempty(file) || ~startsWith(lower(file),lower(root)), continue, end
    end
    % Strip whole-line comments only, so a '%' inside a format string can
    % never hide a later call on the same line (errs toward flagging).
    code=regexprep(fileread(file),'^\s*%[^\n]*','','lineanchors');
    if ~isempty(regexp(code,blocking,'once'))
        offenders{end+1}=name; %#ok<AGROW>
    end
    pending=[pending regexp(code,callPattern,'match')]; %#ok<AGROW>
end
end
