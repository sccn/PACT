function report = validate_pact_filters(varargin)
%VALIDATE_PACT_FILTERS Measure the frequency responses of PACT's FIR paths.
%   report = validate_pact_filters                 % tables to the console
%   report = validate_pact_filters('plots', dir)   % also write PNG figures
%
% For each filter path PACT uses, the filter is applied to a unit impulse
% through the production code (pop_eegfiltnew, pop_firws, filt_fir1fft), and
% the measured impulse response is compared with the design coefficients and
% analysed with FREQZ. Reported per band:
%   alignment  max |measured impulse response - design b| (zero-phase check)
%   symmetric  max |b - flip(b)| (linear phase; 0 means zero group delay
%              after the order/2 compensation)
%   req        requested band edges (Hz)
%   -6 dB      measured -6 dB edges (Hz); -3 dB likewise
%   centre     gain at the geometric band centre (dB)
%   ripple     max |gain| inside the passband (dB); for firws (scanHfo) the
%              passband is [low+df/2, high-df/2]; NaN if the transition fills the band
%   stop       max gain in the stopband beyond the transition band (dB)
%   DC         gain at 0 Hz (dB)
%   sine       amplitude ratio and phase lag (deg) of a centre-frequency
%              sinusoid filtered through the production path, mid-record
%   support    order/2 samples (s): transient length at edges/boundaries
% Paths: (1) pac_man (pop_eegfiltnew); (2) pac_scanHfoByPhase (pop_firws),
% legacy vs revised orders; (3) pac_scanLfoPhaseFreq (filt_fir1fft, Kaiser).
% No random numbers are used.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

plotDir = '';
k = find(strcmpi(varargin,'plots'),1);
if ~isempty(k), plotDir = varargin{k+1}; end
report = struct('path',{},'fs',{},'band',{},'order',{},'alignment',{},'symmetry',{},...
    'minus6',{},'minus3',{},'centreDb',{},'rippleDb',{},'stopDb',{},'dcDb',{},...
    'sineRatio',{},'sineLagDeg',{},'supportSec',{},'warning',{});

fprintf('\n== pac_man filters (pop_eegfiltnew; Hamming; firfilt) ==\n');
configs = {100,[2 4],[20 35]; 200,[3 4],[35 70]; 500,[4 8],[80 150]; 1000,[1 2],[80 150]; 1000,[0.5 1],[100 0]};
for c = 1:size(configs,1)
    fs = configs{c,1};
    for band = configs(c,2:3)
        b = band{1};
        [B,msg] = eegfiltnewCoefficients(fs,b);
        report(end+1) = measure('pac_man',fs,b,B,@(x)eegfiltnewApply(x,fs,b),msg); %#ok<AGROW>
    end
end

fprintf('\n== pac_scanHfoByPhase filters (pop_firws; Hamming) ==\n');
grids = {100,[2 5],[20 35],4; 500,[0.5 8],[100 200],10};
for g = 1:size(grids,1)
    [fs,phaseRange,ampRange,nBins] = grids{g,:};
    for family = {'phase','amp'}
        if strcmp(family{1},'phase')
            list = logspace(log10(phaseRange(1)),log10(phaseRange(2)),nBins); legacyDf = mean(diff(list))/2;
        else
            list = logspace(log10(ampRange(1)),log10(ampRange(2)),nBins); legacyDf = mean(diff(list));
        end
        legacyOrder = ceil(3.3/(legacyDf/fs)/2)*2;
        for i = 1:numel(list)-1
            edges = list(i:i+1);
            for version = {'legacy','revised'}
                if strcmp(version{1},'legacy'), order = legacyOrder; else, order = pac_firwsOrder(fs,edges,legacyDf); end
                [B,msg] = firwsCoefficients(fs,edges,order);
                report(end+1) = measure(['scanHfo/' family{1} '/' version{1}],fs,edges,B,...
                    @(x)firwsApply(x,fs,edges,order),msg); %#ok<AGROW>
            end
        end
    end
end

fprintf('\n== pac_scanLfoPhaseFreq filters (filt_fir1fft; Kaiser, ripple 0.1) ==\n');
fs = 500;
list = logspace(log10(.5),log10(8),15);
bands = [list(1:end-1).' list(2:end).'; 100 0];
for i = 1:size(bands,1)
    edges = bands(i,:);
    [~,order,B] = filt_fir1fft([],edges(1),edges(2),fs);
    report(end+1) = measure('scanLfo/filt_fir1fft',fs,edges,B,...
        @(x)filt_fir1fft(x(:),edges(1),edges(2),fs).',''); %#ok<AGROW>
    assert(numel(B)==order+1);
end

if ~isempty(plotDir)
    plotResponses(report,plotDir);
end
end

% ---------------------------------------------------------------------------
function r = measure(path,fs,band,B,apply,msg)
B = B(:).';
order = numel(B)-1; half = order/2;
n = max(2^nextpow2(4*order+1),2^14);
centre = n/2;
impulse = zeros(1,n); impulse(centre) = 1;
y = evalcQuiet(@()apply(impulse));
measured = y(centre-half:centre+half);
r.path = path; r.fs = fs; r.band = band; r.order = order;
r.alignment = max(abs(measured-B));
r.symmetry = max(abs(B-fliplr(B)));
[H,f] = freqz(B,1,2^18,fs);
gain = 20*log10(abs(H)+realmin);
r.dcDb = gain(1);
lo = band(1); hi = band(2);
if hi == 0, hi = fs/2; end                       % high-pass
if lo == 0, lo = 0; end
r.minus6 = crossings(f,gain,-6.0206,lo,hi);
r.minus3 = crossings(f,gain,-3.0103,lo,hi);
if band(2)==0
    centreF = min(fs/2*.9,2*band(1));
else
    centreF = sqrt(band(1)*band(2));
end
r.centreDb = interp1(f,gain,centreF);
df = 3.3*fs/order;                                 % Hamming nominal width
if startsWith(path,'scanHfo')
    % firws cutoffs are -6 dB points: the passband is the band minus half a
    % transition at each side (empty when the transition fills the band).
    inPass = f>=band(1)+df/2 & f<=band(2)-df/2;
else
    inPass = f>=band(1) & f<=hi;
end
if any(inPass), r.rippleDb = max(abs(gain(inPass))); else, r.rippleDb = NaN; end
if startsWith(path,'scanLfo')                     % Kaiser design: 20% edges
    stop = f<=.8*band(1) | (band(2)>0 & f>=1.2*band(2));
else
    stop = f<=band(1)-df | (band(2)>0 & f>=band(2)+df);
end
r.stopDb = max(gain(stop));
t = (0:max(n,4*order+round(20*fs/centreF))-1)/fs;
x = cos(2*pi*centreF*t);
yx = evalcQuiet(@()apply(x));
mid = round(numel(t)/2)+(-round(fs/centreF):round(fs/centreF));
analyticIn = hilbert(x); analyticOut = hilbert(yx);
r.sineRatio = median(abs(analyticOut(mid))./abs(analyticIn(mid)));
r.sineLagDeg = rad2deg(angle(mean(analyticOut(mid).*conj(analyticIn(mid)))));
r.supportSec = half/fs;
r.warning = msg;
fprintf(['TABLE %-24s fs %5d req [%7.3f %7.3f] order %6d align %.1e sym %.1e ',...
    '-6dB [%7.3f %7.3f] -3dB [%7.3f %7.3f] centre %6.2f dB ripple %5.2f dB stop %7.1f dB ',...
    'DC %7.1f dB sine %.3f/%5.2f deg support %.2f s%s\n'],...
    path,fs,band(1),band(2),order,r.alignment,r.symmetry,r.minus6,r.minus3,r.centreDb,...
    r.rippleDb,r.stopDb,r.dcDb,r.sineRatio,r.sineLagDeg,r.supportSec,ternary(isempty(msg),'',[' WARN: ' msg]));
end

function edges = crossings(f,gain,level,lo,hi)
% Outermost frequencies where the response crosses LEVEL around the band.
above = gain >= level;
centre = find(f>=(lo+hi)/2,1);
if ~above(centre)
    edges = [NaN NaN]; return
end
left = find(~above(1:centre),1,'last');
right = centre-1+find(~above(centre:end),1,'first');
if isempty(left), lowEdge = 0; else, lowEdge = interp1(gain(left:left+1),f(left:left+1),level); end
if isempty(right), highEdge = f(end); else, highEdge = interp1(gain(right-1:right),f(right-1:right),level); end
edges = [lowEdge highEdge];
end

function [B,msg] = eegfiltnewCoefficients(fs,band)
EEG = impulseSet(fs,4096);
evalc('[~,~,B]=pop_eegfiltnew(EEG,band(1),band(2));');
msg = '';
end

function y = eegfiltnewApply(x,fs,band)
EEG = impulseSet(fs,numel(x)); EEG.data = x;
EEG = pop_eegfiltnew(EEG,band(1),band(2));
y = double(EEG.data);
end

function [B,msg] = firwsCoefficients(fs,edges,order)
lastwarn('');
w = warning('off','firfilt:filterOrderLow');
restore = onCleanup(@()warning(w));
B = firws(order,edges/(fs/2),windows('hamming',order+1));
msg = '';
nyquistLimit = min([2*edges 2*(fs/2-edges) diff(edges)]);
if 3.3*fs/order > nyquistLimit
    msg = sprintf('pop_firws: order %d below recommended (transition %.3f Hz > %.3f Hz)',...
        order,3.3*fs/order,nyquistLimit);
end
end

function y = firwsApply(x,fs,edges,order)
EEG = impulseSet(fs,numel(x)); EEG.data = x;
w = warning('off','firfilt:filterOrderLow');
restore = onCleanup(@()warning(w));
EEG = pop_firws(EEG,'fcutoff',edges,'ftype','bandpass','wtype','hamming','forder',order,'minphase',0);
y = double(EEG.data);
end

function EEG = impulseSet(fs,n)
EEG = eeg_emptyset;
EEG.srate = fs; EEG.nbchan = 1; EEG.pnts = n; EEG.trials = 1;
EEG.xmin = 0; EEG.xmax = (n-1)/fs; EEG.data = zeros(1,n);
end

function y = evalcQuiet(fn)
[~,y] = evalc('fn()');
y = double(y(:)).';
end

function out = ternary(condition,a,b)
if condition, out = a; else, out = b; end
end

function plotResponses(report,plotDir)
if ~exist(plotDir,'dir'), mkdir(plotDir); end
visibility = get(groot,'defaultFigureVisible');
restore = onCleanup(@()set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');
groups = {'pac_man','scanHfo/phase','scanHfo/amp','scanLfo'};
for g = 1:numel(groups)
    fig = figure('Position',[0 0 1100 650],'Color','w');
    members = find(startsWith({report.path},groups{g}));
    hold on
    for m = members
        r = report(m);
        if startsWith(r.path,'scanHfo') && r.fs ~= 500, continue, end
        style = '-';
        if contains(r.path,'legacy'), style = ':'; end
        switch true
            case startsWith(r.path,'pac_man'), [B,~] = eegfiltnewCoefficients(r.fs,r.band);
            case startsWith(r.path,'scanHfo'), [B,~] = firwsCoefficients(r.fs,r.band,r.order);
            otherwise, [~,~,B] = filt_fir1fft([],r.band(1),r.band(2),r.fs);
        end
        [H,f] = freqz(B,1,2^16,r.fs);
        plot(f,20*log10(abs(H)+realmin),style,'LineWidth',1);
    end
    yline(-6,'k:'); ylim([-90 5]); set(gca,'XScale','log'); grid on
    xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
    title(sprintf('%s: measured design responses (dotted = legacy order)',groups{g}),'Interpreter','none');
    exportgraphics(fig,fullfile(plotDir,sprintf('filters_%s.png',strrep(groups{g},'/','_'))),'Resolution',110);
    close(fig);
end
end
