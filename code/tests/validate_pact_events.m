function results = validate_pact_events(varargin)
%VALIDATE_PACT_EVENTS Event-mode (hfoPool = 4) trial-shuffle test validation.
%   results = validate_pact_events            % 100 channels x 4 batches
%   results = validate_pact_events('quick')   % 1 batch
%
% End-to-end pac_man runs, fs 200 Hz, LFO 3 Hz (band 2-4 Hz), HFO carrier
% 40-60 Hz (band 35-70 Hz), +/-500 ms windows, HAS 5/20/100 %, 199
% permutations, fixed seeds. Events every 2.5 s (+/-0.3 s jitter) over 120 s
% (about 47 events, Monte Carlo) and, for E07, every 8 s over 60 s (6
% events, exact enumeration over 720 permutations).
%
%   E01 null: events unrelated to the LFO; no coupling             [type I]
%   E02 events at LFO peaks, constant HFO envelope (no modulation)  [type I]
%   E03 identical evoked response at every event (3 Hz ERP + 50 Hz
%       burst at a fixed latency): evoked, not oscillatory PAC      [type I]
%   E04 ongoing PAC inside windows only, k = 0.3; event-random LFO phase [power]
%   E05 as E04 with k = 0.6                                         [power]
%   E06 PAC with phase reset at events (events at LFO peaks), k = 0.6. Phase
%       still varies across events within the window (drift), and the
%       shuffle detects coupling through that variability. Only phase that
%       is identical across events at every latency is undetectable.  [power]
%   E07 E01 with 6 events (exact permutation null)                  [type I]
%   E03c E03 analysed in continuous mode (block sign-flip test)     [contrast]
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

quick = any(strcmpi(varargin,'quick'));
nBatch = 4; if quick, nBatch = 1; end
nCh = 100; fs = 200; rates = [5 20 100]; alpha = .05;
codes = {'E01','E02','E03','E04','E05','E06','E07','E03c'};
rows = {};
for c = 1:numel(codes)
    rej = zeros(1,numel(rates)); tested = rej; missing = rej; minP = nan(1,numel(rates));
    for batch = 1:nBatch
        rng(5000+100*c+batch,'twister');
        EEG = simulate(strtok(codes{c},'c'),fs,nCh);
        for r = 1:numel(rates)
            if strcmp(codes{c},'E03c')
                opts = {'hfoPool',1};
            else
                opts = {'hfoPool',4,'whichMarker',1,'windowLength',500};
            end
            evalc(['o=pac_man(EEG,''lfoPhase'',[2 4],''hfoAmp'',[35 70],''hfoTopRatio'',rates(r),',...
                   '''numSurro'',199,''numPhaseBin'',12,opts{:});']);
            p = o.pac.uncorrected.MIpval;
            ok = isfinite(p);
            rej(r) = rej(r)+sum(p(ok)<=alpha); tested(r) = tested(r)+sum(ok);
            missing(r) = missing(r)+sum(~ok);
            minP(r) = min([minP(r); o.pac.minAttainableP(ok)]);
        end
    end
    for r = 1:numel(rates)
        [rate,ci] = binofit(rej(r),max(tested(r),1));
        rows(end+1,:) = {codes{c},rates(r),rej(r),tested(r),missing(r),100*rate,100*ci(1),100*ci(2),minP(r)}; %#ok<AGROW>
        fprintf('EVENT %-4s HAS%3d%%  rej %4d/%-4d (%5.1f%% [%4.1f,%4.1f])  missing %d  min attainable p %.4g\n',...
            codes{c},rates(r),rej(r),tested(r),100*rate,100*ci(1),100*ci(2),missing(r),minP(r));
    end
end
results = cell2table(rows,'VariableNames',{'scenario','HAS','rejected','tested','missing',...
    'ratePct','ciLowPct','ciHighPct','minAttainableP'});
end

function EEG = simulate(code,fs,nCh)
dur = 120; spacing = 2.5;
if strcmp(code,'E07'), dur = 60; spacing = 8; end
n = dur*fs; t = (0:n-1)/fs;
EEG = eeg_emptyset;
EEG.srate = fs; EEG.pnts = n; EEG.nbchan = nCh; EEG.trials = 1; EEG.xmax = t(end);
shared = any(strcmp(code,{'E02','E06'}));
if shared
    theta = repmat(2*pi*3*t + cumsum(sqrt(.25/fs)*randn(1,n)),nCh,1);
else
    theta = bsxfun(@plus,2*pi*3*t,2*pi*rand(nCh,1)) + cumsum(sqrt(.25/fs)*randn(nCh,n),2);
end
% Event latencies.
if shared
    peaks = find(diff(sign(diff(cos(theta(1,:)))))<0)+1;
    centers = [];
    for target = spacing*fs:spacing*fs:n-spacing*fs
        [~,k] = min(abs(peaks-target-randi(round([-.3 .3]*fs))));
        centers(end+1) = peaks(k); %#ok<AGROW>
    end
    centers = unique(centers);
else
    centers = round(spacing*fs:spacing*fs:n-spacing*fs);
    centers = centers + randi(round([-.3 .3]*fs),size(centers));
end
window = round(.5*fs);
gate = false(1,n);
for c = centers, gate(c-window:c+window) = true; end
% Signals.
phi = exp(-1/(.25*fs));
ar = filter(sqrt(1-phi^2),[1 -phi],randn(nCh,n+500),[],2);
envelope = exp(.5*ar(:,501:end));
spectrum = fft(randn(nCh,n),[],2);
f = abs([0:n/2 -(n/2-1):-1])*fs/n;
carrier = real(ifft(bsxfun(@times,spectrum,f>=40 & f<=60),[],2));
carrier = bsxfun(@rdivide,carrier,std(carrier,0,2));
k = 0;
switch code
    case 'E02'
        envelope(:) = 1;
        carrier = cos(bsxfun(@plus,2*pi*50*t,2*pi*rand(nCh,1)));
    case 'E04', k = .3;
    case 'E05', k = .6;
    case 'E06', k = .6;
end
modulation = 1 + k*bsxfun(@times,gate,cos(theta));
x = cos(theta) + envelope.*modulation.*carrier + .3*randn(nCh,n);
if strcmp(code,'E03')
    tt = (0:round(.8*fs))/fs;
    erp = 1.5*sin(2*pi*3*tt).*exp(-((tt-.25).^2)/(2*.12^2));
    burst = 2*sin(2*pi*50*tt).*exp(-((tt-.15).^2)/(2*.03^2));
    impulses = zeros(1,n); impulses(centers) = 1;
    response = conv(impulses,erp+burst);
    x = x + repmat(response(1:n),nCh,1);
end
EEG.data = x;
EEG.event = struct('type','stim','latency',num2cell(centers),'duration',0);
EEG.urevent = EEG.event;
end
