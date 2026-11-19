function repro_original_cancel_tag_hang
%REPRO_ORIGINAL_CANCEL_TAG_HANG Deliberately reproduce the 18 Nov 2026 hang.
% DO NOT run this interactively or outside run_matlab_with_timeout.ps1: it is
% designed never to return. It reproduces the original harness exactly: a
% timer polls for a uicontrol tagged 'cancel' (case-sensitive) and asserts
% EEGLAB's callback string, while MATLAB's inputdlg (tag 'Cancel') waits in
% uiwait. test_process_timeout.ps1 runs it to show that the process-level
% guard terminates such a hang and reports a timeout instead of waiting.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

dismiss=timer('ExecutionMode','fixedSpacing','Period',.25,'TimerFcn',@originalDismiss);
start(dismiss);
fprintf('REPRO waiting in inputdlg with the original exact-tag harness\n');
inputdlg({'Channel number'},'Event-related PACT',[1 30],{'1'});
fprintf('REPRO unexpectedly returned\n');
end

function originalDismiss(timerObject,~)
% Verbatim logic of the original harness's dismissDialog.
button=findall(groot,'Type','uicontrol','Tag','cancel');
if isempty(button), return, end
assert(strcmp(get(button(1),'Callback'),'close(gcbf)'));
stop(timerObject);
close(ancestor(button(1),'figure'));
end
