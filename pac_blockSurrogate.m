function [stat, nullStat, info] = pac_blockSurrogate(phase, amplitude, selectedIdx, validMask, candidateMask, nSurro, srate, minPhaseHz)
%PAC_BLOCKSURROGATE Block-debiased HAS modulation index and sign-flip null.
%   [stat, nullStat, info] = pac_blockSurrogate(phase, selectedAmplitude, ...
%       HASindices, validMask, candidateMask, nSurro, samplingRate, ...
%       minimumPhaseFrequency)
%
% Inputs
%   phase          1 x T LFO phase (radians) in original sample coordinates.
%   amplitude      1 x N HFO envelope at the HAS samples.
%   selectedIdx    1 x N unique HAS sample indices (subset of candidateMask).
%   validMask      1 x T samples free of edge/boundary filter artifacts.
%   candidateMask  1 x T eligible samples before HAS (validMask restricted to
%                  event windows in event mode).
%   nSurro         number of Monte Carlo sign-flip draws.
%
% Outputs
%   stat      block-debiased modulation index (dMI), defined below.
%   nullStat  1 x nSurro values of dMI under the sign-flip null.
%   info      struct: nBlocks (HAS-occupied blocks), blockLength (samples),
%             blockSums (complex column vector of block contributions).
%
% Method
%   Every uninterrupted run of valid samples is split into non-overlapping
%   blocks of at least max(0.5 s, 5 cycles of minimumPhaseFrequency); a run
%   shorter than that is one block. Blocks never cross EEGLAB boundaries.
%   Within block j, m_j is the mean LFO unit phasor over the candidate
%   samples of that block. Each HAS sample contributes
%       a_t * (exp(i*phase_t) - m_j),
%   and dMI = |sum over all HAS samples| / N.
%
%   Subtracting m_j removes two artifacts of the raw MI that are unrelated
%   to coupling: (1) phase-clustering bias from non-uniform LFO phase
%   distributions, e.g. non-sinusoidal waveforms (cf. van Driel et al.,
%   2015, J Neurosci Methods 254:60-72), and (2) partial-cycle remainders of
%   the mean amplitude, which make adjacent blocks negatively correlated.
%
%   Null hypothesis: the block contributions B_j are independent and
%   centrally symmetric (B_j and -B_j equally likely). This holds
%   approximately when HFO amplitude is independent of LFO phase and blocks
%   are long relative to the dependence of both signals. The null draws
%   independent random signs for the B_j (a block sign-flip randomization
%   test). It does not assume isotropy of B_j. With few occupied blocks the
%   attainable p-values are coarse and the test is conservative.
%
%   This is NOT a trial-level or causal test. Event-locked waveforms that
%   produce both phase concentration and amplitude changes, broadband
%   transients that are filtered into both bands, and slow dependence that
%   crosses block edges violate the null and can yield small p-values.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

phase = phase(:)';
amplitude = amplitude(:)';
selectedIdx = selectedIdx(:)';
validMask = logical(validMask(:)');
candidateMask = logical(candidateMask(:)') & validMask;
nSamples = numel(phase);
if numel(validMask)~=nSamples || numel(candidateMask)~=nSamples || ...
        numel(amplitude)~=numel(selectedIdx)
    error('PACT:InvalidInput','Phase, masks, HAS index and amplitude lengths conflict.');
end
if isempty(selectedIdx) || any(selectedIdx<1 | selectedIdx>nSamples) || ...
        any(selectedIdx~=round(selectedIdx)) || numel(unique(selectedIdx))~=numel(selectedIdx)
    error('PACT:InvalidInput','HAS index must contain unique valid sample indices.');
end
if any(~candidateMask(selectedIdx)) || any(~isfinite(amplitude)) || ...
        any(~isfinite(phase(candidateMask)))
    error('PACT:InvalidInput','HAS includes ineligible samples or nonfinite data.');
end
if ~isscalar(nSurro) || nSurro<1 || nSurro~=round(nSurro) || ...
        ~isscalar(srate) || ~isfinite(srate) || srate<=0 || ...
        ~isscalar(minPhaseHz) || ~isfinite(minPhaseHz) || minPhaseHz<=0
    error('PACT:InvalidInput','Invalid surrogate settings.');
end

blockLength = max(round(0.5*srate),ceil(5*srate/minPhaseHz));
edges = diff([false validMask false]);
segmentStarts = find(edges==1);
segmentEnds = find(edges==-1)-1;
blockId = zeros(1,nSamples);
nAll = 0;
for s = 1:numel(segmentStarts)
    span = segmentStarts(s):segmentEnds(s);
    count = max(1,floor(numel(span)/blockLength));
    % Spread any remainder evenly so that every block has >= blockLength
    % samples unless the whole run is shorter than one block.
    blockId(span) = nAll+1+floor((0:numel(span)-1)*count/numel(span));
    nAll = nAll+count;
end

phasor = exp(1i*phase);
candidates = find(candidateMask);
candidateIds = blockId(candidates);
localMean = accumarray(candidateIds(:),phasor(candidates).',[nAll 1]) ./ ...
            max(1,accumarray(candidateIds(:),1,[nAll 1]));
ids = blockId(selectedIdx);
contribution = amplitude .* (phasor(selectedIdx) - localMean(ids).');
blockSums = accumarray(ids(:),contribution(:),[nAll 1]);
occupied = accumarray(ids(:),1,[nAll 1])>0;
blockSums = blockSums(occupied);
nBlocks = numel(blockSums);
nHAS = numel(selectedIdx);
stat = abs(sum(blockSums))/nHAS;
info = struct('nBlocks',nBlocks,'blockLength',blockLength,'blockSums',blockSums);
if nBlocks<4
    error('PACT:InsufficientBlocks', ...
        'Only %d HAS-containing blocks; at least four are required for exploratory surrogate inference.',nBlocks);
end

% Chunk the draws to cap peak memory; use the caller's RNG state.
nullStat = zeros(1,nSurro);
chunkSize = 256;
for first = 1:chunkSize:nSurro
    last = min(first+chunkSize-1,nSurro);
    signs = 2*(rand(nBlocks,last-first+1)<.5)-1;
    nullStat(first:last) = abs(sum(bsxfun(@times,blockSums,signs),1))/nHAS;
end
end
