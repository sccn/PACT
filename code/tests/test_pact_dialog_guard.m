function test_pact_dialog_guard
%TEST_PACT_DIALOG_GUARD Static check that menu callbacks use only inputgui.
% Regression for the 18 Nov 2026 batch hang: "Plot event-related PAC" called
% MATLAB's inputdlg, whose Cancel button is tagged 'Cancel'. The harness
% searched for 'cancel', never found it, and uiwait blocked indefinitely on an
% invisible figure. The mocked menu test (test_pact_menus) can only answer
% EEGLAB's inputgui, so every dialog reachable from a registered PACT
% callback must use inputgui. This test walks each callback and every PACT
% function it reaches (recursively, PACT files only) and fails on inputdlg,
% questdlg, listdlg, uiwait or keyboard. No figures are drawn.
% Real dialogs are exercised separately in code/tests/gui.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

pactRoot=fileparts(which('eegplugin_pac'));
assert(~isempty(pactRoot),'PACT must be on the MATLAB path.');
callbacks=registeredCallbacks(fullfile(pactRoot,'eegplugin_pac.m'));
assert(numel(callbacks)==11,'Expected 11 registered PACT menu callbacks.');
for k=1:numel(callbacks)
    offenders=pact_findBlockingDialogs(callbacks{k},pactRoot);
    assert(isempty(offenders),'Menu callback "%s" reaches a blocking dialog in: %s',...
        callbacks{k},strjoin(offenders,', '));
end

% The guard must detect the defect it is meant to catch: the 11faa31
% version of pac_vis_eventRelatedPac called inputdlg.
fixture=[tempname '_guard'];
mkdir(fixture);
removeFixture=onCleanup(@()rmdir(fixture,'s'));
fid=fopen(fullfile(fixture,'pac_vis_eventRelatedPac.m'),'w');
fprintf(fid,['function pac_vis_eventRelatedPac(EEG)\n',...
    'answer=inputdlg({''Channel number''},''Event-related PACT'',[1 30],{''1''});\nend\n']);
fclose(fid);
offenders=pact_findBlockingDialogs('pac_vis_eventRelatedPac(EEG);',fixture);
assert(isequal(offenders,{'pac_vis_eventRelatedPac'}),'Guard failed to flag inputdlg.');
fprintf('test_pact_dialog_guard: PASS\n');
end

function callbacks=registeredCallbacks(pluginFile)
code=fileread(pluginFile);
lines=regexp(code,'[^\n]*','match');
callbacks={};
for k=1:numel(lines)
    line=strtrim(lines{k});
    if startsWith(line,'%') || ~contains(line,'uimenu( submenu'), continue, end
    match=regexp(line,'''callback'',\s*''((?:[^'']|'''')*)''','tokens','once');
    if ~isempty(match), callbacks{end+1}=strrep(match{1},'''''',''''); end %#ok<AGROW>
end
end
