function test_pact_integration
%TEST_PACT_INTEGRATION PACT core regression with installed EEGLAB.
% Requires EEGLAB and Signal Processing / Statistics functionality on path.
% No plots, dialogue boxes, or input files are needed.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

assert(exist('eeg_emptyset','file')==2,'EEGLAB must be on MATLAB path.');
pac_setPath;
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
assert(all(isfinite(out.pac.dmi)) && all(out.pac.numBlocks>=4));
assert(isequal(size(out.pac.dmi),[2 1]) && out.pac.blockLength==250);
outRepeat=pac_man(EEG,options{:});
assert(isequaln(out.pac.mi,outRepeat.pac.mi));
assert(isequaln(out.pac.uncorrected.MIpval,outRepeat.pac.uncorrected.MIpval));
for ch=1:2
    ix=out.pac.hfoIndex{ch};
    assert(all(out.pac.validSampleMask(ix)));
    z=out.pac.analyticEEG(ch,ix);
    assert(abs(out.pac.mi(ch)-abs(mean(z)))<1e-12);
    assert(abs(out.pac.angleStd(ch)-circ_std(angle(z(:))))<1e-12);
end

% Event-only analysis excludes all samples outside valid event windows.
options{8}=4; % hfoPool value
options{10}=2; % 'stim' follows numeric boundary -99 in the marker menu.
options{12}=500; % windowLength (ms)
options{6}=100;
centers=500:500:5500;
centers(centers==3000)=[];
centers=centers+round((0:numel(centers)-1)*3.3);
EEG.event=struct('type',{-99,'stim','stim'},'latency',{3000.5,80,3030},'duration',0);
t=(0:n-1)/fs;
lf=cos(2*pi*3*t);
gate=false(1,n);
for center=centers
    EEG.event(end+1)=struct('type','stim','latency',center,'duration',0);
    gate(center-50:center+50)=true;
end
[~,order]=sort([EEG.event.latency]);
EEG.event=EEG.event(order);
EEG.urevent=EEG.event;
EEG.data=repmat(lf+(1+.8*lf.*gate).*cos(2*pi*28*t),2,1);
outEvent=pac_man(EEG,options{:});
assert(isequal(outEvent.data,EEG.data) && isequal(outEvent.event,EEG.event));
assert(isequal(outEvent.urevent,EEG.urevent));
assert(isequal(outEvent.pac.winEdgesAfterRejection,[centers'-50 centers'+50]));
assert(isequal(outEvent.pac.hfoIndex{1},find(gate)));
assert(all(outEvent.pac.uncorrected.MIpval<=.05));
assert(isequal(size(outEvent.pac.windowMeanAmp),[2 101]));
% Event windows use the same angular deviation as angleStd and legacy circ_std.
zWin=exp(1i*angle(outEvent.pac.analyticEEG(1,outEvent.pac.winEdgesAfterRejection(:,1)+50)));
assert(abs(outEvent.pac.windowStdPhase(1,51)-circ_std(angle(zWin(:))))<1e-12);
fprintf('test_pact_integration: PASS\n');
end
