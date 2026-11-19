function [miNull, rNull, ampNull] = pac_shiftSurrogate(phase, amplitude, selectedIdx, validMask, nSurro, srate, nBins)
%PAC_SHIFTSURROGATE Time-shift phase relative to selected amplitudes.
%   Phase is shifted separately within each uninterrupted valid segment;
%   no shifted sample crosses an EEGLAB boundary. The HAS sample indices
%   and observed amplitude values are held fixed. This preserves the
%   temporal autocorrelation of the phase series within each segment,
%   unlike randomly permuting individual phase samples.
%
%   Null quantities: Canolty MI, resultant length, and the range of
%   phase-binned mean amplitudes. Not a replacement for trial-level
%   inference in event-related or nonstationary recordings.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

phase = phase(:)';
amplitude = amplitude(:)';
selectedIdx = selectedIdx(:)';
validMask = logical(validMask(:)');
if numel(phase) ~= numel(validMask)
    error('phase and validMask must have identical lengths.');
end
if any(~validMask(selectedIdx))
    error('HAS samples include invalid/boundary samples.');
end
if any(~isfinite(phase(validMask))) || any(~isfinite(amplitude))
    error('NaN or Inf in phase or selected amplitude.');
end
if nSurro < 1 || nSurro ~= round(nSurro)
    error('nSurro must be a positive integer.');
end

n = numel(selectedIdx);
if n == 0
    error('No HAS samples selected.');
end
d = diff([false validMask false]);
segmentStart = find(d == 1);
segmentEnd = find(d == -1) - 1;

% Precompute per-segment positions and amplitude weights.
groups = {};
for s = 1:numel(segmentStart)
    sel = find(selectedIdx >= segmentStart(s) & selectedIdx <= segmentEnd(s));
    if isempty(sel)
        continue
    end
    segmentLength = segmentEnd(s) - segmentStart(s) + 1;
    if segmentLength < 4
        error('Selected samples in a segment shorter than four samples.');
    end
    entry.phase = phase(segmentStart(s):segmentEnd(s));
    entry.positions = selectedIdx(sel) - segmentStart(s) + 1;
    entry.outputIdx = sel;
    entry.amplitude = amplitude(sel);
    entry.length = segmentLength;
    % Keep at least one sample of separation, and prefer >= 1 s
    % where the uninterrupted segment is long enough.
    entry.minShift = min(max(1, round(srate)), max(1, floor(segmentLength/4)));
    groups{end+1} = entry; %#ok<AGROW>
end

miNull = zeros(1, nSurro);
rNull = zeros(1, nSurro);
ampNull = nan(1, nSurro);
edges = linspace(-pi, pi, nBins+1);
for iteration = 1:nSurro
    phaseSamples = zeros(1,n);
    for s = 1:numel(groups)
        g = groups{s};
        shift = randi([g.minShift, g.length-g.minShift]);
        shifted = mod(g.positions-1+shift, g.length) + 1;
        phaseSamples(g.outputIdx) = g.phase(shifted);
    end
    z = exp(1i*phaseSamples);
    miNull(iteration) = abs(mean(amplitude .* z));
    rNull(iteration) = abs(mean(z));
    bin = discretize(phaseSamples, edges);
    binMean = accumarray(bin(:), amplitude(:), [nBins 1], @mean, NaN);
    if sum(isfinite(binMean)) >= 2
        ampNull(iteration) = max(binMean) - min(binMean);
    end
end
end
