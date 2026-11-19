function test_pact_scans
%TEST_PACT_SCANS Independent recomputation of both PACT frequency scans.
% Requires EEGLAB (firfilt). Deterministic (fixed seed). No dialogs; figures
% are drawn invisibly where the scan always plots.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

visibility = get(groot,'defaultFigureVisible');
restore = onCleanup(@()set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');
existing = findall(groot,'Type','figure');
closeNew = onCleanup(@()delete(setdiff(findall(groot,'Type','figure'),existing)));

rng(61,'twister');
fs = 200; n = 120*fs; t = (0:n-1)/fs;
EEG = eeg_emptyset;
EEG.srate = fs; EEG.nbchan = 2; EEG.pnts = n; EEG.trials = 1; EEG.xmax = t(end);
theta = 2*pi*3*t + cumsum(sqrt(.25/fs)*randn(1,n));
EEG.data = [cos(theta)+(1+.5*cos(theta)).*cos(2*pi*45*t)+.3*randn(1,n);
            randn(1,n)];
% Mixed marker types: a numeric -99 boundary and a text event.
EEG.event = struct('type',{'stim',-99},'latency',{1000,12000.5},'duration',0);
EEG.urevent = EEG.event;
original = EEG;

% ---------------- pac_scanHfoByPhase ----------------
phaseRange = [2 6]; ampRange = [30 60]; nBins = 4; HAS = 5;
evalc('out = pac_scanHfoByPhase(EEG,phaseRange,ampRange,HAS,nBins,3,1);');
assert(isequal(out.data,original.data) && isequaln(out.event,original.event),'Scan changed EEG.');
S = out.pacScan;
assert(isequal(size(S.lfoPhaseMi),[2 nBins-1 nBins-1]) && isequal(size(S.resulVectLen),size(S.lfoPhaseMi)));
phaseList = logspace(log10(phaseRange(1)),log10(phaseRange(2)),nBins);
ampList = logspace(log10(ampRange(1)),log10(ampRange(2)),nBins);
assert(max(abs(S.phaseFreqEdges-phaseList))<1e-12 && max(abs(S.ampFreqEdges-ampList))<1e-12);
% Every band's transition fits inside the band (no collapsed passbands).
for k = 1:nBins-1
    assert(3.3*fs/S.phaseFilterOrder(k) <= diff(phaseList(k:k+1))+1e-9);
    assert(3.3*fs/S.ampFilterOrder(k) <= diff(ampList(k:k+1))+1e-9);
end
% Independent recomputation: firfilt-style constant padding per segment,
% segment-wise Hilbert, guard covering half the longest filter.
cuts = [1 12001 n+1];                          % boundary at 12000.5 -> 12001
guard = max(fs,max([S.phaseFilterOrder S.ampFilterOrder])/2);
valid = false(1,n);
for s = 1:2
    valid(cuts(s)+guard:cuts(s+1)-1-guard) = true;
end
assert(isequal(valid,S.validSampleMask) && S.guardSamples==guard);
for pair = [1 1; 2 3; 3 2].'
    pk = pair(1); ak = pair(2);
    phase = analytic(EEG.data,cuts,valid,phaseList(pk:pk+1),S.phaseFilterOrder(pk),fs,'angle');
    amp = analytic(EEG.data,cuts,valid,ampList(ak:ak+1),S.ampFilterOrder(ak),fs,'abs');
    for ch = 1:2
        a = amp(ch,:); p = phase(ch,:);
        sel = a > prctile(a,100-HAS);
        mi = abs(mean(a(sel).*exp(1i*p(sel))));
        rvl = abs(mean(exp(1i*p(sel))));
        assert(abs(double(S.lfoPhaseMi(ch,pk,ak))-mi) <= 1e-5*max(mi,1e-3),'MI cell mismatch.');
        assert(abs(double(S.resulVectLen(ch,pk,ak))-rvl) <= 1e-5,'RVL cell mismatch.');
        assert(abs(double(S.meanHfoAmp(ch,pk,ak))-mean(a(sel))) <= 1e-5*mean(a(sel)));
    end
end
% The coupled channel shows the expected maximum at the 3 Hz / 45 Hz cell.
[~,best] = max(reshape(S.lfoPhaseMi(1,:,:),1,[]));
[bp,ba] = ind2sub([nBins-1 nBins-1],best);
assert(phaseList(bp)<=3 && phaseList(bp+1)>=3 && ampList(ba)<=45 && ampList(ba+1)>=45);
% 100% HAS uses every valid sample.
evalc('all100 = pac_scanHfoByPhase(EEG,phaseRange,ampRange,100,nBins,3,1);');
amp = analytic(EEG.data,cuts,valid,ampList(1:2),all100.pacScan.ampFilterOrder(1),fs,'abs');
assert(abs(double(all100.pacScan.meanHfoAmp(2,1,1))-mean(amp(2,:))) <= 1e-5*mean(amp(2,:)));
% Unsupported: filters longer than the data segments are rejected.
short = EEG; short.pnts = 4*fs; short.data = EEG.data(:,1:short.pnts); short.xmax = 4-1/fs;
short.event = []; short.urevent = [];
expectError(@()quiet(@()pac_scanHfoByPhase(short,[0.5 1],ampRange,HAS,4,3,1)),'PACT:ScanUnsupported');

% ---------------- pac_scanLfoPhaseFreq ----------------
noBoundary = EEG; noBoundary.event = noBoundary.event(1); noBoundary.urevent = noBoundary.event;
% Regression: a numeric -99 boundary used to slip past the boundary check.
expectError(@()quiet(@()pac_scanLfoPhaseFreq(EEG,[2 6],4,[30 60],[1 10],3,1,1)),'');
evalc('lfo = pac_scanLfoPhaseFreq(noBoundary,[2 6],4,[30 60],[1 10],3,1,1);');
assert(isequal(lfo.data,original.data),'LFO scan changed EEG.data.');
L = lfo.pacScan;
assert(isequal(size(L.lfoPhaseMi),[2 3 3]));
guard = L.guardSamples; validIdx = guard+1:n-guard;
[~,~,bAmp] = filt_fir1fft([],30,60,fs);
amp = abs(hilbert(reflectFilter(double(EEG.data).',bAmp))).';
amp = amp(:,validIdx);
rates = logspace(log10(1),log10(10),3);
phaseEdges = logspace(log10(2),log10(6),4);
for pk = [1 3]
    [~,~,bPhase] = filt_fir1fft([],phaseEdges(pk),phaseEdges(pk+1),fs);
    phase = angle(hilbert(reflectFilter(double(EEG.data).',bPhase))).';
    phase = phase(:,validIdx);
    for r = [1 3]
        count = ceil(rates(r)/100*numel(validIdx));
        for ch = 1:2
            [~,order] = sort(amp(ch,:),'descend');
            top = order(1:count);
            mi = abs(mean(amp(ch,top).*exp(1i*phase(ch,top))));
            assert(abs(L.lfoPhaseMi(ch,pk,r)-mi) <= 1e-9*max(mi,1),'LFO scan cell mismatch.');
        end
    end
end

% ---------------- runtime (informational) ----------------
bench = EEG; bench.event = []; bench.urevent = [];
bench.data = repmat(EEG.data,4,1); bench.nbchan = 8;
% Timer handles: firfilt calls a bare TIC internally, which resets TOC.
clock = tic; evalc('pac_scanHfoByPhase(bench,[0.5 8],[30 90],2,10,3,1);'); hfoSeconds = toc(clock);
clock = tic; evalc('pac_scanLfoPhaseFreq(bench,[0.5 8],15,[30 80],[0.1 3],15,1,1);'); lfoSeconds = toc(clock);
expectError(@()quiet(@()pac_scanLfoPhaseFreq(bench,[0.5 8],15,[30 90],[0.1 3],15,1,1)),'PACT:ScanUnsupported');
fprintf('SCAN BENCH 8 ch x 120 s @ 200 Hz: HFOxPhase 10x10 grid %.1f s; HASxLFO 15x15 grid %.1f s\n',...
    hfoSeconds,lfoSeconds);
fprintf('test_pact_scans: PASS\n');
end

function out = analytic(data,cuts,valid,edges,order,fs,part)
% firfilt arithmetic: pad each segment with order/2 copies of its first and
% last sample, convolve, keep the centred part; Hilbert within the segment.
b = firws(order,edges/(fs/2),windows('hamming',order+1));
half = order/2;
out = zeros(size(data,1),nnz(valid));
column = 0;
for s = 1:numel(cuts)-1
    span = cuts(s):cuts(s+1)-1;
    keep = valid(span);
    if ~any(keep), continue, end
    for ch = 1:size(data,1)
        x = double(data(ch,span));
        y = conv([repmat(x(1),1,half) x repmat(x(end),1,half)],b,'valid');
        h = hilbert(y);
        if strcmp(part,'abs'), v = abs(h); else, v = angle(h); end
        out(ch,column+(1:nnz(keep))) = v(keep);
    end
    column = column+nnz(keep);
end
end

function y = reflectFilter(x,b)
% filt_fir1fft arithmetic: point reflection about each endpoint, zero phase.
n = numel(b)-1;
prefix = 2*x(1,:) - x(n+1:-1:2,:);
suffix = 2*x(end,:) - x(end-1:-1:end-n,:);
padded = [prefix; x; suffix];
full = zeros(size(padded,1)+n,size(x,2));
for c = 1:size(x,2), full(:,c) = conv(padded(:,c),b(:)); end
y = full(n+n/2+1:n+n/2+size(x,1),:);
end

function quiet(fn)
evalc('fn();');
end

function expectError(callback,identifier)
try
    callback();
catch ME
    assert(isempty(identifier) || strcmp(ME.identifier,identifier),'Unexpected error %s: %s',ME.identifier,ME.message);
    return
end
error('PACT:Test','Expected an error did not occur.');
end
