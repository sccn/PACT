function test_pact_integration
%TEST_PACT_INTEGRATION PACT core regression with installed EEGLAB.
% Requires EEGLAB and Signal Processing / Statistics functionality on path.
% No plots, dialogue boxes, or input files are needed.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

assert(exist('eeg_emptyset','file')==2,'EEGLAB must be on MATLAB path.');
rng(4);
fs=100;
n=6000;
EEG=eeg_emptyset;
EEG.srate=fs;
EEG.data=randn(2,n);
EEG.pnts=n;
EEG.nbchan=2;
EEG.trials=1;
EEG.xmin=0;
EEG.xmax=(n-1)/fs;
EEG.event=struct('type',{'stim','boundary','stim'},...
                 'latency',{2000,3000,4500},'duration',{0,0,0});
EEG.urevent=EEG.event;
EEG.chanlocs=struct('labels',{'Ch1','Ch2'});
inputData=EEG.data;
inputEvents=EEG.event;
options={'lfoPhase',[2 4],'hfoAmp',[20 35],...
    'hfoTopRatio',5,'hfoPool',1,'whichMarker',1,...
    'windowLength',[],'alpha',.05,'numSurro',50,'numPhaseBin',12};
out=pac_man(EEG,options{:});
assert(isequal(out.data,inputData));
assert(isequal(out.event,inputEvents));
assert(isequal(size(out.pac.analyticEEG),size(inputData)));
assert(~out.pac.validSampleMask(3000));
assert(all(isfinite(out.pac.mi)));
assert(all(out.pac.uncorrected.MIpval>0));
outRepeat=pac_man(EEG,options{:});
assert(isequaln(out.pac.mi,outRepeat.pac.mi));
assert(isequaln(out.pac.uncorrected.MIpval,outRepeat.pac.uncorrected.MIpval));
for ch=1:2
    assert(all(out.pac.validSampleMask(out.pac.hfoIndex{ch})));
end

% Event-only analysis excludes all samples outside valid event windows.
options{8}=4; % hfoPool value
options{12}=100; % windowLength (ms)
outEvent=pac_man(EEG,options{:});
assert(isequal(outEvent.data,inputData));
assert(size(outEvent.pac.winEdgesAfterRejection,1)==2);
assert(all(outEvent.pac.selectedSampleMask(outEvent.pac.hfoIndex{1})));
assert(isequal(size(outEvent.pac.windowMeanAmp),[2 21]));
fprintf('test_pact_integration: PASS\n');
end
