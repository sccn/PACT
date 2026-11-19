function smoke_pact_real_dialogs
%SMOKE_PACT_REAL_DIALOGS Cancel every real PACT menu dialog (no mocks).
% Real-GUI smoke test, kept separate from the mocked callback test
% (code/tests/test_pact_menus.m). Each registered PACT menu that opens a
% dialog is invoked through its actual callback string with EEGLAB's real
% inputgui. pact_dialogDismisser presses the real Cancel button. The test
% fails if a dialog is not reached, is not canceled before the per-dialog
% deadline, or if EEG/ALLEEG change. Run it in its own MATLAB process under
% code/tests/gui/run_matlab_with_timeout.ps1, which adds an independent
% wall-clock limit for failures that MATLAB itself cannot recover from.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

deadline=20;
existing=findall(groot,'Type','figure');
closeNew=onCleanup(@()delete(setdiff(findall(groot,'Type','figure'),existing)));
host=figure('Visible','off');
uimenu(host,'tag','tools');
eegplugin_pac(host,struct,struct);
menu=findall(host,'type','uimenu','label','PACT');
items=flipud(get(menu,'Children'));
assert(numel(items)==11);

rng(79);
EEG=eeg_emptyset; EEG.srate=100; EEG.nbchan=2; EEG.pnts=8000; EEG.trials=1;
EEG.xmax=79.99; EEG.data=randn(2,8000);
EEG.chanlocs=struct('labels',{'A','B'});
EEG.event=struct('type','stim','latency',num2cell(500:500:7500),'duration',0);
EEG.urevent=EEG.event;
% Compute PAC non-interactively so that statistics and event dialogs exist.
evalc(['EEG=pac_man(EEG,''lfoPhase'',[2 4],''hfoAmp'',[20 35],''hfoTopRatio'',20,',...
       '''hfoPool'',4,''whichMarker'',1,''windowLength'',500,''numSurro'',49,''numPhaseBin'',12);']);
CURRENTSET=1; ALLEEG=[EEG EEG]; %#ok<NASGU>

withDialog={'Compute PAC','Plot HFO-marked raw data','Copy event markers',...
    'Set up statistics','Plot event-related PAC','Scan HAS x Phase_Freq',...
    'Scan HFO_Freq x Phase_Freq'};
labels=arrayfun(@(h)get(h,'label'),items,'UniformOutput',false);
for k=1:numel(withDialog)
    item=items(strcmp(labels,withDialog{k}));
    assert(isscalar(item),'Menu "%s" not registered.',withDialog{k});
    before=EEG; beforeAll=ALLEEG;
    dismisser=pact_dialogDismisser(deadline);
    started=tic;
    try
        eval(get(item,'callback'));
    catch ME
        pact_dialogDismisserStop(dismisser);
        rethrow(ME);
    end
    elapsed=toc(started);
    log=pact_dialogDismisserStop(dismisser);
    assert(~isempty(log),'Menu "%s" did not open a dialog.',withDialog{k});
    actions={log.action};
    assert(~any(strcmp(actions,'timeout')),'Menu "%s": dialog timeout (%s).',...
        withDialog{k},log(find(strcmp(actions,'timeout'),1)).detail);
    assert(isequaln(before,EEG) && isequaln(beforeAll,ALLEEG),...
        'Menu "%s": canceling changed EEG or ALLEEG.',withDialog{k});
    fprintf('REAL CANCEL PASS %-28s %.2f s\n',withDialog{k},elapsed);
end
fprintf('smoke_pact_real_dialogs: PASS\n');
end
