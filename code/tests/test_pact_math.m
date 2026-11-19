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
% Surrogate draws are reproducible for a fixed caller RNG.
fs=100;
t=(1:1000)/fs;
phase=2*pi*2*t;
valid=true(1,1000);
valid(450:550)=false;
selected=[200:210 800:810];
amps=1+0.1*cos(phase(selected));
rng(9);
a=pac_shiftSurrogate(phase,amps,selected,valid,50,fs,12);
rng(9);
b=pac_shiftSurrogate(phase,amps,selected,valid,50,fs,12);
assert(isequal(a,b));
assert(numel(a)==50 && all(isfinite(a)));
assert(all(a>=0));
fprintf('test_pact_math: PASS\n');
end
