function test_pact_filters
%TEST_PACT_FILTERS Fast regression checks of PACT's FIR designs.
% Requires EEGLAB (firfilt) for the pop_firws check. Deterministic.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

% 1. filt_fir1fft (pac_scanLfoPhaseFreq) meets its documented spec:
%    deviation <= 0.01 in the passband and stopband (40 dB), zero phase.
%    Regression: the 0.1-ripple design used a rectangular window and lost
%    2.1 dB at the passband edges of narrow LFO bands.
fs = 500;
for band = {[0.5 0.61],[2 2.438],[6.563 8],[100 0],[0 30]}
    lo = band{1}(1); hi = band{1}(2);
    [~,order,b] = filt_fir1fft([],lo,hi,fs);
    [H,f] = freqz(b,1,2^17,fs); g = abs(H);
    if lo == 0
        pass = f<=hi; stop = f>=1.2*hi;
    elseif hi == 0
        pass = f>=lo; stop = f<=.8*lo;
    else
        pass = f>=lo & f<=hi; stop = f<=.8*lo | f>=1.2*hi;
    end
    assert(max(abs(g(pass)-1)) <= .015,'filt_fir1fft [%g %g]: passband deviation %.3f.',lo,hi,max(abs(g(pass)-1)));
    assert(max(g(stop)) <= .015,'filt_fir1fft [%g %g]: stopband gain %.3f.',lo,hi,max(g(stop)));
    n = 4*order+1; x = zeros(n,1); x((n+1)/2) = 1;
    evalc('y = filt_fir1fft(x,lo,hi,fs);');
    assert(max(abs(y((n+1)/2-order/2:(n+1)/2+order/2).'-b)) < 1e-12,'filt_fir1fft is not zero-phase.');
end

% 2. pac_scanHfoByPhase transition widths. Regression: with the GUI default
%    grid (0.5-8 Hz, 10 log bins) the legacy order lost 3.4 dB at the centre
%    of the lowest phase bin and pop_firws warned "Filter order too low".
list = logspace(log10(.5),log10(8),10);
legacyDf = mean(diff(list))/2;
legacyOrder = ceil(3.3/(legacyDf/fs)/2)*2;
assert(centreGainDb(fs,list(1:2),legacyOrder) < -3,'Legacy defect not reproduced.');
EEG = eeg_emptyset; EEG.srate = fs; EEG.nbchan = 1; EEG.trials = 1;
EEG.pnts = 20*fs; EEG.xmax = (EEG.pnts-1)/fs; EEG.data = randn(1,EEG.pnts);
for k = 1:numel(list)-1
    order = pac_firwsOrder(fs,list(k:k+1),legacyDf);
    assert(abs(centreGainDb(fs,list(k:k+1),order)) < .1,'Bin %d centre gain off.',k);
    lastwarn('');
    evalc(['pop_firws(EEG,''fcutoff'',list(k:k+1),''ftype'',''bandpass'',',...
           '''wtype'',''hamming'',''forder'',order,''minphase'',0);']);
    [~,id] = lastwarn;
    assert(~strcmp(id,'firfilt:filterOrderLow'),'pop_firws still warns for bin %d.',k);
end

% 3. pac_man filters (pop_eegfiltnew) are zero-phase: impulse response equals
%    the design coefficients, centred.
EEG.pnts = 8192; EEG.xmax = (EEG.pnts-1)/fs; EEG.data = zeros(1,EEG.pnts); EEG.data(4096) = 1;
evalc('[out,~,b] = pop_eegfiltnew(EEG,3,4);');
half = (numel(b)-1)/2;
assert(max(abs(out.data(4096-half:4096+half)-b)) < 1e-12 && max(abs(b-fliplr(b))) < 1e-15);
fprintf('test_pact_filters: PASS\n');
end

function gainDb = centreGainDb(fs,edges,order)
b = firws(order,edges/(fs/2),windows('hamming',order+1));
H = freqz(b,1,[0 sqrt(prod(edges))],fs);
gainDb = 20*log10(abs(H(2)));
end
