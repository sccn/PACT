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

% Regression (Nov 2026): non-sinusoidal LFO with independent HFO amplitude.
% The block-rotation null rejected ~100% of such channels at 20% HAS.
rng(31,'twister');
nCh=40; n=9000; t=(0:n-1)/fs;
null=eeg_emptyset; null.srate=fs; null.pnts=n; null.nbchan=nCh; null.trials=1;
null.xmax=t(end);
theta=bsxfun(@plus,2*pi*3*t,2*pi*rand(nCh,1))+cumsum(.03*randn(nCh,n),2);
ar=filter(sqrt(1-.96^2),[1 -.96],randn(nCh,n+2000),[],2);
null.data=cos(theta)+.6*cos(2*theta+pi/3)+.3*cos(3*theta+1)+...
    exp(.45*ar(:,2001:end)).*cos(bsxfun(@plus,2*pi*28*t,2*pi*rand(nCh,1)))+.15*randn(nCh,n);
evalc(['outNull=pac_man(null,''lfoPhase'',[2 10],''hfoAmp'',[20 36],',...
    '''hfoTopRatio'',20,''numSurro'',199);']);
assert(sum(outNull.pac.uncorrected.MIpval<=.05)<=6,...
    'Non-sinusoidal null rejected %d of 40 channels.',sum(outNull.pac.uncorrected.MIpval<=.05));

% Regression: events locked to one LFO phase with a constant HFO envelope.
% Raw MI is large, but there is no amplitude modulation; the rotation null
% gave p = 0.005. Block-local debiasing must not reject.
rng(3,'twister');
fs2=200; n2=fs2*120; t2=(0:n2-1)/fs2;
locked=eeg_emptyset; locked.srate=fs2; locked.pnts=n2; locked.nbchan=1;
locked.trials=1; locked.xmax=t2(end);
theta=2*pi*4*t2+2*pi*rand;
locked.data=cos(theta)+cos(2*pi*28*t2+2*pi*rand)+.05*randn(1,n2);
peaks=find(diff(sign(diff(cos(theta))))<0)+1;
peaks=peaks(peaks>2*fs2 & peaks<n2-2*fs2);
locked.event=struct('type','stim','latency',num2cell(peaks(1:8:end)),'duration',0);
locked.urevent=locked.event;
evalc(['outLocked=pac_man(locked,''lfoPhase'',[3 5],''hfoAmp'',[20 36],''hfoTopRatio'',100,',...
    '''hfoPool'',4,''whichMarker'',1,''windowLength'',40,''numSurro'',199,''numPhaseBin'',4);']);
assert(outLocked.pac.mi>.5 && outLocked.pac.uncorrected.MIpval>.05);
% Regression: 25-s blocks (75 cycles of 3 Hz) leave 4 blocks in 120 s (p >= 0.125).
% The p-value is valid but can never reach 0.05; pac_man must warn.
lastwarn('');
evalc(['short=pac_man(locked,''lfoPhase'',[3 5],''hfoAmp'',[20 36],''hfoTopRatio'',20,',...
    '''numSurro'',99,''blockCycles'',75,''numPhaseBin'',4);']);
[~,id]=lastwarn;
assert(short.pac.exactNull && short.pac.minAttainableP>.05 && strcmp(id,'PACT:Underpowered'));
fprintf('test_pact_integration: PASS\n');
end
