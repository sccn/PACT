function test_pact_math
%TEST_PACT_MATH Unit tests; no EEGLAB, GUI or EEG datasets required.
% MATLAB command:
%   addpath('/path/to/PACT'); test_pact_math
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

p=[.001 .01 .02 .2 NaN];
holm=pac_adjustPValues(p,'holm');
expected=[.004 .03 .04 .2 NaN];
assert(max(abs(holm(1:4)-expected(1:4))) < 1e-12);
assert(isnan(holm(5)));
bh=pac_adjustPValues(p,'fdr');
expectedBH=[.004 .02 .02*4/3 .2 NaN];
assert(max(abs(bh(1:4)-expectedBH(1:4))) < 1e-12);
assert(isequal(pac_adjustPValues([0 1],'bonferroni'),[0 1]));
assert(isequal(pac_parseNumeric('[1 2.5e-3]'),[1 .0025]));
assert(isempty(pac_parseNumeric('')));
try
    pac_parseNumeric('system(''echo bad'')');
    error('PACT:Test','Expected parsing failure did not occur.');
catch ME
    assert(~strcmp(ME.identifier,'PACT:Test'));
end
% Block-debiased statistic and sign-flip null (pac_blockSurrogate).
fs=100;
t=(1:3000)/fs;
phase=2*pi*2*t+.3*sin(2*pi*.37*t);
valid=true(1,3000);
valid(1450:1550)=false;
selected=[200:210 450:460 700:710 900:910 1800:1810 2100:2110 2400:2410 2700:2710];
amps=1+0.1*cos(phase(selected));
rng(9);
[s1,a,info]=pac_blockSurrogate(phase,amps,selected,valid,valid,50,fs,2);
rng(9);
[s2,b]=pac_blockSurrogate(phase,amps,selected,valid,valid,50,fs,2);
assert(isequal(a,b) && s1==s2);
% Eight occupied blocks: exact enumeration of all 2^7 sign patterns.
assert(info.exact && numel(a)==128 && all(isfinite(a)) && all(a>=0));
assert(info.blockLength==250 && info.nBlocks==8);
assert(abs(info.minAttainableP-2^-7)<1e-15 && info.pval>=2^-7 && info.pval<=1);
assert(abs(info.pval-mean(a>=s1-1e-12))<1e-12);
% Independent reference: two valid runs (1449 and 1450 samples), five
% blocks each, remainder spread evenly; subtract each block's mean phasor.
ids=zeros(1,3000);
ids(1:1449)=1+floor((0:1448)*5/1449);
ids(1551:3000)=6+floor((0:1449)*5/1450);
z=exp(1i*phase);
reference=0;
for blk=unique(ids(selected))
    m=mean(z(ids==blk));
    hit=ids(selected)==blk;
    reference=reference+sum(amps(hit).*(z(selected(hit))-m));
end
assert(abs(s1-abs(reference)/numel(selected))<1e-12);
assert(all(a<=sum(abs(info.blockSums))/numel(selected)+1e-12));
% A global phase offset cannot change the statistic.
assert(abs(pac_blockSurrogate(phase+.7,amps,selected,valid,valid,5,fs,2)-s1)<1e-12);
% Constant amplitude at 100% selection is exactly unbiased: dMI = 0,
% even for a non-uniform (here strongly skewed) phase distribution.
skewed=phase+.8*sin(phase);
allValid=find(valid);
[zeroStat,~,zeroInfo]=pac_blockSurrogate(skewed,ones(size(allValid)),allValid,valid,valid,5,fs,2);
assert(zeroStat==0 && zeroInfo.pval==1,'Degenerate zero statistic must give p = 1.');
% blockCycles sets the block length: 10 cycles of 2 Hz at 100 Hz = 500 samples.
[~,~,info10]=pac_blockSurrogate(phase,amps,selected,valid,valid,5,fs,2,10);
assert(info10.blockLength==500);
% More than 16 occupied blocks: Monte Carlo with (1+b)/(N+1).
wide=1:60:3000; wide=wide(valid(wide));
[~,mc,infoMC]=pac_blockSurrogate(phase,1+cos(phase(wide)),wide,valid,valid,99,fs,20,.5);
assert(~infoMC.exact && infoMC.nBlocks>16 && numel(mc)==99 && infoMC.minAttainableP==1/100);
% Runs shorter than one block form their own block instead of failing.
short=false(1,3000);
short([1:120 300:420 600:720 900:1020])=true;
[~,~,info]=pac_blockSurrogate(phase,ones(1,4),[60 360 660 960],short,short,5,fs,2);
assert(info.nBlocks==4);
expectError(@()pac_blockSurrogate(phase,ones(1,3),[60 360 660],short,short,5,fs,2),...
    'PACT:InsufficientBlocks');
candidate=valid;
candidate(200)=false;
expectError(@()pac_blockSurrogate(phase,amps,selected,valid,candidate,5,fs,2),...
    'PACT:InvalidInput');
% Event shuffle (pac_eventShuffle): exact for W <= 8, independent reference.
rng(5);
T=2000; ph=2*pi*rand(1,T); am=1+rand(1,T);
win=[(101:200:1501)' (101:200:1501)'+49];            % 8 windows of 50 samples
sel=sort(reshape(bsxfun(@plus,win(:,1),[0 10 20 30 40]).',1,[]));
[es,enull,einfo]=pac_eventShuffle(ph,am,sel,win,10);
assert(einfo.exact && numel(enull)==factorial(8) && abs(einfo.minAttainableP-1/factorial(8))<1e-15);
assert(abs(es-abs(mean(am(sel).*exp(1i*ph(sel)))))<1e-12);
order=perms(1:8); pick=order(1234,:); ref=0;
for w=1:8
    offs=sel(sel>=win(w,1) & sel<=win(w,2))-win(w,1);
    ref=ref+sum(am(win(w,1)+offs).*exp(1i*ph(win(pick(w),1)+offs)));
end
assert(any(abs(enull*numel(sel)-abs(ref))<1e-9),'Shuffle null misses an enumerated permutation.');
expectError(@()pac_eventShuffle(ph,am,sel,[win(1:7,:); win(7,:)+10],10),'PACT:OverlappingWindows');
expectError(@()pac_eventShuffle(ph,am,sel(1:15),win(1:3,:),10),'PACT:InsufficientEvents');
fprintf('test_pact_math: PASS\n');
end

function expectError(callback,identifier)
try
    callback();
catch ME
    assert(strcmp(ME.identifier,identifier),'Unexpected error %s.',ME.identifier);
    return
end
error('PACT:Test','Expected error %s did not occur.',identifier);
end
