function test_pact_dialog_timeout
%TEST_PACT_DIALOG_TIMEOUT Real-dialog regression for the Cancel-tag hang.
% On 18 Nov 2026 a batch menu test waited indefinitely: a timer searched for
% a uicontrol tagged 'cancel', but MATLAB's inputdlg tags its button 'Cancel'
% (findall is case-sensitive), so uiwait never returned. This test opens
% real MATLAB dialogs (no PACT callbacks, no mocks) and shows that:
%   1. the original exact-tag search cannot find inputdlg's button;
%   2. pact_dialogDismisser cancels the same inputdlg through its own
%      callback, so inputdlg returns {} promptly;
%   3. a modal dialog with no recognizable Cancel button is deleted at the
%      deadline and reported as 'timeout', so no wait can be indefinite.
% Run under code/tests/gui/run_matlab_with_timeout.ps1 for an independent
% process-level limit.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

existing=findall(groot,'Type','figure');
cleanupFigures=onCleanup(@()delete(setdiff(findall(groot,'Type','figure'),existing)));
% 1. Root cause, observed on a live inputdlg: an inspection timer records
%    which tags exist, then deletes the figure so this step always ends.
setappdata(groot,'PACTtagProbe',[]);
probe=timer('ExecutionMode','fixedSpacing','Period',.2,'TasksToExecute',100,...
    'TimerFcn',@inspectInputdlg);
removeProbe=onCleanup(@()delete(probe));
start(probe);
started=tic;
answer=inputdlg({'Channel number'},'Event-related PACT',[1 30],{'1'});
assert(toc(started)<20,'Root-cause probe did not end promptly.');
tags=getappdata(groot,'PACTtagProbe');
assert(~isempty(tags),'inputdlg was not reached by the probe.');
assert(tags.lowercaseFound==false,'Original harness tag search unexpectedly matched.');
assert(tags.capitalizedFound==true,'inputdlg no longer uses the Cancel tag.');
assert(isempty(answer));
fprintf('ROOT CAUSE reproduced: exact tag ''cancel'' absent, ''Cancel'' present.\n');

% 2. The revised harness cancels inputdlg via the button's own callback.
dismisser=pact_dialogDismisser(10);
started=tic;
answer=inputdlg({'Channel number'},'Event-related PACT',[1 30],{'1'});
elapsed=toc(started);
log=pact_dialogDismisserStop(dismisser);
assert(isempty(answer),'inputdlg did not return a cancellation.');
assert(numel(log)==1 && strcmp(log(1).action,'canceled'),'Expected one cancellation.');
assert(elapsed<5,'inputdlg cancellation took %.1f s.',elapsed);
fprintf('inputdlg canceled by revised harness in %.2f s.\n',elapsed);

% 3. A modal dialog without any Cancel button is bounded by the deadline.
deadline=2;
dismisser=pact_dialogDismisser(deadline);
started=tic;
f=figure('Name','No cancel button','WindowStyle','normal');
uicontrol(f,'Style','pushbutton','String','Abort','Tag','abortButton');
uiwait(f);
elapsed=toc(started);
log=pact_dialogDismisserStop(dismisser);
assert(~isvalid(f),'Dialog was not removed.');
assert(numel(log)==1 && strcmp(log(1).action,'timeout'),'Timeout was not reported.');
assert(elapsed>=deadline && elapsed<deadline+5,'Deadline not honored (%.1f s).',elapsed);
fprintf('Unanswerable dialog removed after %.2f s and reported as timeout.\n',elapsed);
clear cleanupFigures removeProbe
fprintf('test_pact_dialog_timeout: PASS\n');
end

function inspectInputdlg(timerObject,~)
f=findall(groot,'Type','figure','Name','Event-related PACT');
% Inspect only once inputdlg is blocked in uiwait (fully built).
if isempty(f) || ~strcmp(get(f(1),'WaitStatus'),'waiting'), return, end
stop(timerObject);
tags.lowercaseFound=~isempty(findall(f(1),'Type','uicontrol','Tag','cancel'));
tags.capitalizedFound=~isempty(findall(f(1),'Type','uicontrol','Tag','Cancel'));
setappdata(groot,'PACTtagProbe',tags);
delete(f(1));
end
