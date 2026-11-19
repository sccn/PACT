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
assert(numel(a)==50 && all(isfinite(a)) && all(a>=0));
assert(info.blockLength==250 && info.nBlocks==8);
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
assert(pac_blockSurrogate(skewed,ones(size(allValid)),allValid,valid,valid,5,fs,2)<1e-12);
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
