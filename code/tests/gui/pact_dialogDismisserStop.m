function log = pact_dialogDismisserStop(dismisser)
%PACT_DIALOGDISMISSERSTOP Disarm a pact_dialogDismisser timer; return its log.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.
log = struct('name',{},'action',{},'seconds',{},'detail',{});
if isempty(dismisser) || ~isvalid(dismisser), return, end
stop(dismisser);
log = dismisser.UserData.log;
delete(dismisser);
end
