function [stat, nullStat, info] = pac_eventShuffle(phase, amplitude, selectedIdx, windows, nSurro)
%PAC_EVENTSHUFFLE Event-mode HAS MI with a trial (window) shuffle null.
%   [stat, nullStat, info] = pac_eventShuffle(phase, amplitude, HASindices, ...
%       windows, nSurro)
%
% Inputs
%   phase, amplitude  1 x T LFO phase (radians) and HFO envelope of one
%                     channel in original sample coordinates.
%   selectedIdx       HAS sample indices; all must lie inside the windows.
%   windows           W x 2 inclusive [first last] sample of each event
%                     window; all windows have the same length L and must
%                     not overlap.
%   nSurro            Monte Carlo permutations (used when W > 8).
%
% Statistic: the HAS modulation index restricted to the windows,
%   stat = |sum_w sum_l H(w,l) a(w,l) exp(i phase(w,l))| / N,
% where H marks HAS samples and N is their number (equal to pac_man's mi).
%
% Null: the phase of window w is replaced by the phase of window pi(w) at the
% same latency, for random permutations pi; amplitudes and HAS positions
% stay fixed. This tests whether HFO amplitude covaries with LFO phase
% across events beyond what is common to all events. Event-locked (evoked)
% waveforms that are the same in every event, including phase concentration
% with an unmodulated envelope and evoked HFO bursts at a fixed evoked
% phase, are preserved by the shuffle and are therefore not detected as
% coupling. Conversely, genuine PAC with complete phase reset at every event
% is indistinguishable from an evoked response and is not detected either.
% Exchangeability of events requires non-overlapping windows.
%
% P-value: exact over all W! permutations when W <= 8, otherwise
% (1 + #{null >= stat})/(nSurro + 1). info.minAttainableP gives the smallest
% possible p (1/W! when exact). At least four events are required.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

phase = phase(:).';
amplitude = amplitude(:).';
selectedIdx = selectedIdx(:).';
if size(windows,2)~=2 || isempty(windows) || any(windows(:)~=round(windows(:))) || ...
        any(windows(:,2)<windows(:,1)) || any(windows(:)<1) || any(windows(:)>numel(phase))
    error('PACT:InvalidInput','windows must be W x 2 inclusive sample ranges.');
end
lengths = windows(:,2)-windows(:,1)+1;
if any(lengths~=lengths(1))
    error('PACT:InvalidInput','All event windows must have the same length.');
end
nWin = size(windows,1);
if nWin < 4
    error('PACT:InsufficientEvents','Only %d event windows; at least four are required.',nWin);
end
ordered = sortrows(windows);
if any(ordered(2:end,1) <= ordered(1:end-1,2))
    error('PACT:OverlappingWindows', ...
        'Event windows overlap; the trial-shuffle null needs non-overlapping windows.');
end
if ~isscalar(nSurro) || nSurro<1 || nSurro~=round(nSurro)
    error('PACT:InvalidInput','Invalid number of permutations.');
end

L = lengths(1);
index = bsxfun(@plus,windows(:,1),0:L-1);           % W x L sample indices
isHAS = false(1,numel(phase));
isHAS(selectedIdx) = true;
if nnz(isHAS(index)) ~= numel(selectedIdx)
    error('PACT:InvalidInput','Every HAS sample must lie inside an event window.');
end
weighted = amplitude(index).*isHAS(index);          % W x L
phasor = exp(1i*phase(index));                      % W x L
pairs = weighted*phasor.';                          % pairs(w,v): amp of w, phase of v
nHAS = numel(selectedIdx);
stat = abs(trace(pairs))/nHAS;
tolerance = 1e-10*max(sum(abs(pairs(:))),eps);
info = struct('nEvents',nWin,'pval',NaN,'exact',false,'minAttainableP',NaN);

if nWin <= 8
    order = perms(1:nWin);
    linear = sub2ind([nWin nWin],repmat(1:nWin,size(order,1),1),order);
    nullStat = abs(sum(pairs(linear),2)).'/nHAS;
    info.pval = mean(nullStat*nHAS >= stat*nHAS - tolerance);
    info.exact = true;
    info.minAttainableP = 1/size(order,1);
    return
end
nullStat = zeros(1,nSurro);
for s = 1:nSurro
    nullStat(s) = abs(sum(pairs(sub2ind([nWin nWin],1:nWin,randperm(nWin)))))/nHAS;
end
info.pval = (1+sum(nullStat*nHAS >= stat*nHAS - tolerance))/(nSurro+1);
info.minAttainableP = 1/(nSurro+1);
end
