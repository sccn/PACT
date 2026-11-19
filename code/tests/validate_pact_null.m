function results = validate_pact_null(varargin)
%VALIDATE_PACT_NULL Adversarial validation of the continuous-mode dMI test.
%   results = validate_pact_null                    % development grid
%   results = validate_pact_null('holdout')         % holdout grid
%   results = validate_pact_null('centering')       % phase-mean variants
%   results = validate_pact_null(..., 'quick')      % 1 batch, smoke run
%   results = validate_pact_null(..., 'csv', file)  % also write a CSV
%   results = validate_pact_null(..., 'only', {'N02'}, 'batches', 16, 'seedOffset', 50000)
%                                                   % targeted replication
%
% Independent of the production surrogate code in three ways:
%   * own signal generator (scenarios below; fixed seeds);
%   * own reference implementation of the block-debiased statistic and the
%     sign-flip test (refTest, written without calling PACT functions),
%     run on batch 1 and compared with production pac_blockSurrogate;
%   * production results come from pac_man (filtering, masks) followed by
%     pac_blockSurrogate for each HAS rate and block length.
%
% Grid (development): HAS 1/2/5/20/100 % x blockCycles 3/5/10/20 x
% recording length 60/300 s, 100 channels x 4 batches per cell, fs 200 Hz,
% LFO 3 Hz (phase band 2-4 Hz), HFO carrier 40-60 Hz (amplitude band 35-70).
% Holdout: fs 500 Hz, LFO 6 Hz (band 4-8), carrier 90-130 Hz (band 80-150),
% 120 s, blockCycles 5 and 10, different seeds; scenarios as listed below.
%
% Scenario codes (k = coupling strength)
%   N01 stationary sinusoidal LFO, AR(1) log-envelope (tau 0.25 s)
%   N02 non-sinusoidal LFO (harmonics), phase band 2-10 Hz
%   N03 nonstationary LFO frequency (2.5-3.5 Hz) and amplitude
%   N04 pink (1/f) noise only
%   N05 long-memory (1/f) log-envelope
%   N06 on/off HFO bursts with gamma-distributed durations
%   N07 strongly skewed envelope (log-normal, sigma 1)
%   N08 slow gain drift and a gain step (whole signal)
%   N09 six independently generated segments joined at EEGLAB boundaries
%   N10 intermittent 50 Hz narrowband artifact bursts
%   N11 broadband single-sample impulses of random sign
%   N12 unequal block occupancy: four 10x HFO bursts per channel
%   C01 sharp-edged sawtooth LFO (no envelope coupling)          [dependence]
%   C02 common-input transients: spike-and-wave shape at random times [dependence]
%   C03 event-locked evoked LFO wave + HFO burst (continuous mode) [dependence]
%   P01 coupling k = 0.10;  P02 k = 0.30
%   P03 k = 0.30 with preferred phase rotating once per minute
% N scenarios have independent LFO phase and HFO amplitude: rejections are
% type-I errors. C scenarios contain genuine statistical dependence that is
% not oscillatory physiological PAC; rejections are expected and show what
% the test cannot exclude. P scenarios measure power.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

args = varargin;
if any(strcmpi(args,'centering'))
    results = centeringStudy(any(strcmpi(args,'quick')));
    return
end
holdout = any(strcmpi(args,'holdout'));
quick = any(strcmpi(args,'quick'));
csvFile = '';
k = find(strcmpi(args,'csv'),1);
if ~isempty(k), csvFile = args{k+1}; end

if holdout
    cfg = struct('fs',500,'f0',6,'lfo',[4 8],'lfoWide',[4 16],'carrier',[90 130],...
        'hfo',[80 150],'durations',120,'cycles',[5 10],'seedBase',9000);
    codes = {'N01','N02','N05','N06','N09','N12','P01'};
else
    cfg = struct('fs',200,'f0',3,'lfo',[2 4],'lfoWide',[2 10],'carrier',[40 60],...
        'hfo',[35 70],'durations',[60 300],'cycles',[3 5 10 20],'seedBase',1000);
    codes = {'N01','N02','N03','N04','N05','N06','N07','N08','N09','N10','N11',...
        'N12','C01','C02','C03','P01','P02','P03'};
end
rates = [1 2 5 20 100];
nCh = 100; nBatch = 4; if quick, nBatch = 1; end
k = find(strcmpi(args,'only'),1);
if ~isempty(k), codes = args{k+1}; end
k = find(strcmpi(args,'batches'),1);
if ~isempty(k), nBatch = args{k+1}; end
k = find(strcmpi(args,'seedOffset'),1);
if ~isempty(k), cfg.seedBase = cfg.seedBase + args{k+1}; end
nSurro = 199; alpha = .05; familySize = 10; minHAS = 36;

rows = {};
agreement = struct('statMaxDiff',0,'exactMismatch',0,'exactCompared',0,'mcAbsDiff',[]);
dependence = {};
for c = 1:numel(codes)
    for d = 1:numel(cfg.durations)
        dur = cfg.durations(d);
        tally = zeros(numel(rates),numel(cfg.cycles),7); % rej,tested,missing,famRej,fams,refRej,refTested
        small = zeros(numel(rates),numel(cfg.cycles),3);  % K in 4..8: tested, rejected; minP>alpha
        adjCos = [];
        for batch = 1:nBatch
            seed = cfg.seedBase + 100*c + 10*d + batch;
            rng(seed,'twister');
            EEG = simulate(codes{c},cfg,dur,nCh);
            band = cfg.lfo;
            if strcmp(codes{c},'N02'), band = cfg.lfoWide; end
            evalc(['o=pac_man(EEG,''lfoPhase'',band,''hfoAmp'',cfg.hfo,''hfoTopRatio'',100,',...
                   '''numSurro'',1,''numPhaseBin'',4);']);
            z = o.pac.analyticEEG; valid = o.pac.validSampleMask; cand = o.pac.selectedSampleMask;
            candIdx = find(cand);
            if batch==1 && dur==cfg.durations(1)
                checkEndToEnd(EEG,band,cfg,o);
            end
            for r = 1:numel(rates)
                for b = 1:numel(cfg.cycles)
                    p = nan(nCh,1); pRef = nan(nCh,1);
                    for ch = 1:nCh
                        a = abs(z(ch,:)); ph = angle(z(ch,:));
                        if rates(r)==100
                            sel = candIdx;
                        else
                            sel = candIdx(a(candIdx) > prctile(a(candIdx),100-rates(r)));
                        end
                        if numel(sel) < minHAS, continue, end
                        rng(seed+ch,'twister');
                        try
                            [st,~,info] = pac_blockSurrogate(ph,a(sel),sel,valid,cand,nSurro,...
                                cfg.fs,band(1),cfg.cycles(b));
                            p(ch) = info.pval;
                            if info.nBlocks <= 8
                                small(r,b,1) = small(r,b,1)+1;
                                small(r,b,2) = small(r,b,2)+(info.pval<=alpha);
                                small(r,b,3) = small(r,b,3)+(info.minAttainableP>alpha);
                            end
                        catch ME
                            if ~strcmp(ME.identifier,'PACT:InsufficientBlocks'), rethrow(ME), end
                            info = [];
                        end
                        if batch == 1
                            [rs,rp,ri] = refTest(ph,a,sel,valid,cand,cfg.fs,band(1),...
                                cfg.cycles(b),nSurro,seed+ch);
                            pRef(ch) = rp;
                            if ~isempty(info)
                                agreement.statMaxDiff = max(agreement.statMaxDiff,abs(rs-st));
                                if info.exact
                                    agreement.exactCompared = agreement.exactCompared+1;
                                    agreement.exactMismatch = agreement.exactMismatch+(abs(rp-info.pval)>1e-12);
                                else
                                    agreement.mcAbsDiff(end+1) = abs(rp-info.pval);
                                end
                            end
                            if rates(r)==100 && cfg.cycles(b)==5 && ~isempty(ri.adjacentCos)
                                adjCos(end+1) = mean(ri.adjacentCos); %#ok<AGROW>
                            end
                        end
                    end
                    ok = isfinite(p);
                    tally(r,b,1) = tally(r,b,1)+sum(p(ok)<=alpha);
                    tally(r,b,2) = tally(r,b,2)+sum(ok);
                    tally(r,b,3) = tally(r,b,3)+sum(~ok);
                    for f = 1:familySize:nCh
                        adj = pac_adjustPValues(p(f:f+familySize-1),'holm');
                        if any(isfinite(adj))
                            tally(r,b,5) = tally(r,b,5)+1;
                            tally(r,b,4) = tally(r,b,4)+any(adj<=alpha);
                        end
                    end
                    okr = isfinite(pRef);
                    tally(r,b,6) = tally(r,b,6)+sum(pRef(okr)<=alpha);
                    tally(r,b,7) = tally(r,b,7)+sum(okr);
                end
            end
        end
        for r = 1:numel(rates)
            for b = 1:numel(cfg.cycles)
                t = squeeze(tally(r,b,:));
                [rate,ci] = binofit(t(1),max(t(2),1));
                [refRate,refCi] = binofit(t(6),max(t(7),1));
                rows(end+1,:) = {codes{c},dur,rates(r),cfg.cycles(b),t(1),t(2),t(3),...
                    100*rate,100*ci(1),100*ci(2),t(4),t(5),100*t(4)/max(t(5),1),...
                    100*refRate,100*refCi(1),100*refCi(2),small(r,b,1),small(r,b,2),small(r,b,3)}; %#ok<AGROW>
                fprintf(['RESULT %s %3ds HAS%3g%% cyc%2d  rej %4d/%-4d (%5.1f%% [%4.1f,%4.1f]) ',...
                    'missing %4d  HolmFWER %3d/%-3d  ref %5.1f%% [%4.1f,%4.1f]  K<=8 %d (rej %d, minP>a %d)\n'],...
                    codes{c},dur,rates(r),cfg.cycles(b),t(1),t(2),100*rate,100*ci(1),100*ci(2),...
                    t(3),t(4),t(5),100*refRate,100*refCi(1),100*refCi(2),small(r,b,:));
            end
        end
        dependence(end+1,:) = {codes{c},dur,mean(adjCos)}; %#ok<AGROW>
        fprintf('ADJACENT %s %3ds HAS100%% cyc5 mean cos(B_j,B_j+1) = %+.3f over %d channels\n',...
            codes{c},dur,mean(adjCos),numel(adjCos));
    end
end
fprintf('AGREEMENT stat max |ref-prod| = %.3g; exact p mismatches %d of %d; MC median |dp| = %.4f (n=%d)\n',...
    agreement.statMaxDiff,agreement.exactMismatch,agreement.exactCompared,...
    median(agreement.mcAbsDiff),numel(agreement.mcAbsDiff));
results = cell2table(rows,'VariableNames',{'scenario','seconds','HAS','blockCycles',...
    'rejected','tested','missing','ratePct','ciLowPct','ciHighPct','familyRejected',...
    'families','fwerPct','refRatePct','refCiLowPct','refCiHighPct','smallK','smallKRejected',...
    'smallKUntestable'});
if ~isempty(csvFile), writetable(results,csvFile); end
end

% ---------------------------------------------------------------------------
function T = centeringStudy(quick)
% Candidate-sample vs HAS-only block phase means vs no centering.
cfg = struct('fs',200,'f0',3,'lfo',[2 4],'lfoWide',[2 10],'carrier',[40 60],...
    'hfo',[35 70],'seedBase',7000);
codes = {'N01','N02','P01','P02'}; rates = [1 5 20 100];
variants = {'candidate','has','none'};
nCh = 100; nBatch = 2; if quick, nBatch = 1; end
rows = {};
for c = 1:numel(codes)
    rej = zeros(numel(rates),numel(variants)); tested = rej;
    for batch = 1:nBatch
        seed = cfg.seedBase+100*c+batch; rng(seed,'twister');
        EEG = simulate(codes{c},cfg,60,nCh);
        band = cfg.lfo; if strcmp(codes{c},'N02'), band = cfg.lfoWide; end
        evalc(['o=pac_man(EEG,''lfoPhase'',band,''hfoAmp'',cfg.hfo,''hfoTopRatio'',100,',...
               '''numSurro'',1,''numPhaseBin'',4);']);
        z = o.pac.analyticEEG; cand = o.pac.selectedSampleMask; candIdx = find(cand);
        for r = 1:numel(rates)
            for ch = 1:nCh
                a = abs(z(ch,:)); ph = angle(z(ch,:));
                if rates(r)==100, sel = candIdx; else, sel = candIdx(a(candIdx) > prctile(a(candIdx),100-rates(r))); end
                for v = 1:numel(variants)
                    [~,p] = refTest(ph,a,sel,o.pac.validSampleMask,cand,cfg.fs,band(1),5,199,seed+ch,variants{v});
                    if isfinite(p), rej(r,v) = rej(r,v)+(p<=.05); tested(r,v) = tested(r,v)+1; end
                end
            end
        end
    end
    for r = 1:numel(rates)
        for v = 1:numel(variants)
            rows(end+1,:) = {codes{c},rates(r),variants{v},rej(r,v),tested(r,v),100*rej(r,v)/max(tested(r,v),1)}; %#ok<AGROW>
        end
        fprintf('CENTERING %s HAS%3d%%  candidate %5.1f%%  HAS-only %5.1f%%  none %5.1f%%  (n=%d)\n',codes{c},rates(r),...
            100*rej(r,1)/max(tested(r,1),1),100*rej(r,2)/max(tested(r,2),1),100*rej(r,3)/max(tested(r,3),1),tested(r,1));
    end
end
T = cell2table(rows,'VariableNames',{'scenario','HAS','centering','rejected','tested','ratePct'});
end

% ---------------------------------------------------------------------------
function checkEndToEnd(EEG,band,cfg,o)
% pac_man at 5% HAS / 5 cycles must give the same dmi as the grid's path.
evalc(['e2e=pac_man(EEG,''lfoPhase'',band,''hfoAmp'',cfg.hfo,''hfoTopRatio'',5,',...
       '''numSurro'',19,''numPhaseBin'',36,''blockCycles'',5);']);
z = o.pac.analyticEEG; cand = find(o.pac.selectedSampleMask);
for ch = 1:min(5,EEG.nbchan)
    a = abs(z(ch,:)); sel = cand(a(cand) > prctile(a(cand),95));
    % |amp*exp(i*phase)| can differ from amp by rounding, which may move a
    % sample lying exactly on the percentile threshold.
    assert(numel(setxor(sel,e2e.pac.hfoIndex{ch}))<=2,'HAS selection differs from pac_man.');
    if isequal(sel,e2e.pac.hfoIndex{ch}) && isfinite(e2e.pac.dmi(ch))
        st = pac_blockSurrogate(angle(z(ch,:)),a(sel),sel,o.pac.validSampleMask,...
            o.pac.selectedSampleMask,19,cfg.fs,band(1),5);
        assert(abs(st-e2e.pac.dmi(ch))<=1e-10*max(1,st),'dmi differs from pac_man.');
    end
end
end

% ---------------------------------------------------------------------------
function [stat,p,info] = refTest(phase,amp,sel,valid,cand,fs,fmin,cycles,nSurro,seed,centering)
% Independent reference: explicit loops, own RNG stream; no PACT calls.
% centering: 'candidate' (production: mean phasor of the block's candidate
% samples), 'has' (mean phasor of the block's HAS samples) or 'none'.
if nargin < 11, centering = 'candidate'; end
info.adjacentCos = [];
L = max(round(.5*fs),ceil(cycles*fs/fmin));
blockOf = zeros(size(valid));
next = 0; i = 1; T = numel(valid);
while i <= T
    if ~valid(i), i = i+1; continue, end
    j = i; while j < T && valid(j+1), j = j+1; end
    runLength = j-i+1; nb = max(1,floor(runLength/L));
    for s = i:j
        blockOf(s) = next + 1 + floor((s-i)*nb/runLength);
    end
    next = next+nb; i = j+1;
end
isSel = false(size(valid)); isSel(sel) = true;
sums = zeros(next,1); occupied = false(next,1);
for blk = 1:next
    inBlock = blockOf==blk;
    c = inBlock & cand;
    if ~any(c), continue, end
    h = inBlock & isSel;
    switch centering
        case 'candidate', m = sum(exp(1i*phase(c)))/nnz(c);
        case 'has',       m = sum(exp(1i*phase(h)))/max(nnz(h),1);
        otherwise,        m = 0;
    end
    if any(h)
        occupied(blk) = true;
        sums(blk) = sum(amp(h).*(exp(1i*phase(h))-m));
    end
end
if all(occupied) && next > 1
    u = sums./max(abs(sums),eps);
    info.adjacentCos = real(u(1:end-1).*conj(u(2:end)));
end
B = sums(occupied); K = numel(B); N = numel(sel);
if sum(abs(B)) <= 1e-10*sum(amp(sel)), B(:) = 0; end
stat = abs(sum(B))/N;
if K < 4, p = NaN; return, end
if ~any(B), p = 1; return, end
tol = 1e-10*sum(abs(B));
if K <= 16
    count = 0; total = 2^(K-1);
    for code = 0:total-1
        s = ones(K,1);
        for bit = 1:K-1
            if bitget(code,bit), s(bit+1) = -1; end
        end
        count = count + (abs(sum(s.*B)) >= abs(sum(B)) - tol);
    end
    p = count/total;
else
    stream = RandStream('mlfg6331_64','Seed',seed);
    exceed = 0;
    for draw = 1:nSurro
        s = sign(rand(stream,K,1)-.5); s(s==0) = 1;
        exceed = exceed + (abs(sum(s.*B)) >= abs(sum(B)) - tol);
    end
    p = (1+exceed)/(nSurro+1);
end
end

% ---------------------------------------------------------------------------
function EEG = simulate(code,cfg,dur,nCh)
fs = cfg.fs; n = round(fs*dur); t = (0:n-1)/fs;
EEG = eeg_emptyset;
EEG.srate = fs; EEG.pnts = n; EEG.nbchan = nCh; EEG.trials = 1; EEG.xmin = 0; EEG.xmax = t(end);
EEG.chanlocs = struct('labels',arrayfun(@(c)sprintf('S%d',c),1:nCh,'UniformOutput',false));
if strcmp(code,'N09')
    seg = 6; len = floor(n/seg); EEG.data = zeros(nCh,n);
    for s = 1:seg
        span = (s-1)*len + (1:len);
        if s == seg, span = (s-1)*len+1:n; end
        EEG.data(:,span) = base(cfg,t(1:numel(span)),nCh,0,.25,.5);
    end
    EEG.event = struct('type','boundary','latency',num2cell((1:seg-1)*len+.5),'duration',0);
    EEG.urevent = EEG.event;
    return
end
k = 0; tau = .25; sigma = .5;
switch code
    case 'N05', tau = -1;                 % long-memory envelope (see base)
    case 'N07', sigma = 1;
    case 'P01', k = .10;
    case 'P02', k = .30;
    case 'P03', k = .30;
end
[x,theta,extra] = base(cfg,t,nCh,k,tau,sigma,code);
switch code
    case 'N04'
        x = pinkNoise(nCh,n,1);
    case 'N06'
        x = x + extra.carrier.*(4*burstTrain(nCh,n,fs,.5,[2 .25]));
    case 'N08'
        gain = 1+.5*sin(bsxfun(@plus,2*pi*t/40,2*pi*rand(nCh,1)));
        stepAt = randi([round(n/4) round(3*n/4)],nCh,1);
        gain(bsxfun(@ge,1:n,stepAt)) = 1.8*gain(bsxfun(@ge,1:n,stepAt));
        x = x.*gain;
    case 'N10'
        on = burstTrain(nCh,n,fs,.1,[3 .6]);
        x = x + 2*on.*sin(bsxfun(@plus,2*pi*50*t,2*pi*rand(nCh,1)));
    case 'N11'
        for c = 1:nCh
            idx = randi([1 n],1,round(8*dur/60));
            x(c,idx) = x(c,idx) + 8*sign(randn(1,numel(idx)));
        end
    case 'N12'
        for c = 1:nCh
            for b = 1:4
                first = randi([1 n-round(.3*fs)]); span = first:first+round(.3*fs)-1;
                x(c,span) = x(c,span) + 9*extra.envelope(c,span).*extra.carrier(c,span);
            end
        end
    case 'C01'
        saw = mod(theta,2*pi)/pi-1;               % sharp falling edge each cycle
        x = x - cos(theta) + 1.5*saw;
    case 'C02'
        x = x + transients(nCh,n,fs,.5);
    case 'C03'
        [x2,events] = evoked(nCh,n,fs);
        x = x + x2;
        EEG.event = events; EEG.urevent = events;
end
EEG.data = x;
end

function [x,theta,extra] = base(cfg,t,nCh,k,tau,sigma,code)
if nargin < 7, code = 'N01'; end
fs = cfg.fs; n = numel(t);
f0 = cfg.f0*ones(1,n);
if strcmp(code,'N03'), f0 = cfg.f0*(1+(1/6)*sin(2*pi*t/50)); end
drift = cumsum(sqrt(.25/fs)*randn(nCh,n),2);
theta = bsxfun(@plus,2*pi*cumsum(f0)/fs,2*pi*rand(nCh,1)) + drift;
lfoAmp = ones(nCh,n);
if strcmp(code,'N03'), lfoAmp = 1+.7*sin(bsxfun(@plus,2*pi*t/30,2*pi*rand(nCh,1))); end
lf = lfoAmp.*cos(theta);
if strcmp(code,'N02'), lf = lf + .6*cos(2*theta+pi/3) + .3*cos(3*theta+1); end
if tau > 0
    phi = exp(-1/(tau*fs));
    ar = filter(sqrt(1-phi^2),[1 -phi],randn(nCh,n+round(10*tau*fs)),[],2);
    logEnv = ar(:,end-n+1:end);
else
    logEnv = pinkNoise(nCh,n,1);
end
envelope = exp(sigma*logEnv);
carrier = bandNoise(nCh,n,fs,cfg.carrier);
preferred = zeros(1,n);
if strcmp(code,'P03'), preferred = 2*pi*t/60; end
modulation = 1 + k*cos(bsxfun(@minus,theta,preferred));
x = lf + envelope.*modulation.*carrier + .3*pinkNoise(nCh,n,1);
extra.carrier = carrier; extra.envelope = envelope;
end

function x = pinkNoise(nCh,n,alpha)
spectrum = fft(randn(nCh,n),[],2);
f = [0:floor(n/2) -(ceil(n/2)-1):-1]; f(1) = 1;
x = real(ifft(bsxfun(@times,spectrum,abs(f).^(-alpha/2)),[],2));
x = bsxfun(@rdivide,bsxfun(@minus,x,mean(x,2)),std(x,0,2));
end

function x = bandNoise(nCh,n,fs,band)
spectrum = fft(randn(nCh,n),[],2);
f = abs([0:floor(n/2) -(ceil(n/2)-1):-1])*fs/n;
x = real(ifft(bsxfun(@times,spectrum,f>=band(1) & f<=band(2)),[],2));
x = bsxfun(@rdivide,x,std(x,0,2));
end

function on = burstTrain(nCh,n,fs,rate,gammaShapeScale)
% On/off bursts: Poisson onsets (per second), gamma durations (seconds).
on = zeros(nCh,n);
for c = 1:nCh
    count = poissrnd(rate*n/fs);
    for b = 1:count
        len = max(1,round(gamrnd(gammaShapeScale(1),gammaShapeScale(2))*fs));
        first = randi([1 n]); on(c,first:min(n,first+len-1)) = 1;
    end
end
kernel = hann(round(.05*fs)).'; kernel = kernel/sum(kernel);
on = conv2(on,kernel,'same');
end

function x = transients(nCh,n,fs,rate)
% Spike-and-wave transients (same polarity) at Poisson times.
tt = (-round(.3*fs):round(.5*fs))/fs;
shape = 4*exp(-tt.^2/(2*.008^2)) - 1.5*exp(-(tt-.1).^2/(2*.08^2));
x = zeros(nCh,n);
for c = 1:nCh
    times = randi([1 n],1,poissrnd(rate*n/fs));
    impulses = zeros(1,n); impulses(times) = 1;
    x(c,:) = conv(impulses,shape,'same');
end
end

function [x,events] = evoked(nCh,n,fs)
% Identical evoked response at every event: 3 Hz Gabor ERP + HFO burst.
period = round(2.5*fs);
centers = period:period:n-period;
centers = centers + randi(round([-.3 .3]*fs),size(centers));
tt = (0:round(.8*fs))/fs;
erp = 1.5*sin(2*pi*3*tt).*exp(-((tt-.25).^2)/(2*.12^2));
burst = 2*sin(2*pi*50*tt).*exp(-((tt-.15).^2)/(2*.03^2));
impulses = zeros(1,n); impulses(centers) = 1;
response = conv(impulses,erp+burst);
x = repmat(response(1:n),nCh,1);
events = struct('type','stim','latency',num2cell(centers),'duration',0);
end
