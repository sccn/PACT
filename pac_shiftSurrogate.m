function miNull = pac_shiftSurrogate(phase, amplitude, selectedIdx, validMask, nSurro, srate, minPhaseHz)
%PAC_SHIFTSURROGATE Block-phase-rotation surrogate for pooled PACT MI.
%   miNull = pac_shiftSurrogate(phase, selectedAmplitude, HASindices, ...
%       validMask, nSurro, samplingRate, minimumPhaseFrequency)
%
% We test *cross-block consistency of preferred coupling phase*, conditional
% on the observed within-block amplitude-phase relationships. Each valid,
% uninterrupted segment is divided into non-overlapping blocks spanning
% at least five cycles at the lowest LFO frequency (minimum 0.5 second).
% Each block's complex PAC sum is independently rotated by a uniformly
% drawn angle, leaving its magnitude and within-block temporal structure
% unchanged. No block crosses an EEGLAB boundary.
%
% This is NOT a generic test of phase-amplitude independence, and does
% not support inference for one-off events or nonstationary transitions.
% At least four HAS-containing blocks are required.
%
% Whole-recording circular shifts are deliberately NOT used: for a pure
% sinusoid, they change only global phase and cannot alter |MI|.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

phase = phase(:)';
amplitude = amplitude(:)';
selectedIdx = selectedIdx(:)';
validMask = logical(validMask(:)');
if numel(phase)~=numel(validMask) || numel(amplitude)~=numel(selectedIdx)
    error('PACT:InvalidInput','Phase, valid mask, HAS index and amplitude lengths conflict.');
end
if isempty(selectedIdx) || any(selectedIdx<1 | selectedIdx>numel(phase)) || ...
        any(selectedIdx~=round(selectedIdx)) || numel(unique(selectedIdx))~=numel(selectedIdx)
    error('PACT:InvalidInput','HAS index must contain unique valid sample indices.');
end
if any(~validMask(selectedIdx)) || any(~isfinite(amplitude)) || ...
        any(~isfinite(phase(validMask)))
    error('PACT:InvalidInput','HAS includes invalid samples or nonfinite data.');
end
if ~isscalar(nSurro) || nSurro<1 || nSurro~=round(nSurro) || ...
        ~isscalar(srate) || ~isfinite(srate) || srate<=0 || ...
        ~isscalar(minPhaseHz) || ~isfinite(minPhaseHz) || minPhaseHz<=0
    error('PACT:InvalidInput','Invalid surrogate settings.');
end

blockLength = max(round(0.5*srate),round(5*srate/minPhaseHz));
segments = diff([false validMask false]);
segmentStarts = find(segments==1);
segmentEnds = find(segments==-1)-1;

% Each block contributes an observed complex sum; block count matters,
% not the number of autocorrelated samples inside a high-amplitude burst.
blockSums = complex(zeros(0,1));
for s=1:numel(segmentStarts)
    first = segmentStarts(s);
    last = segmentEnds(s);
    while first<=last
        final=min(first+blockLength-1,last);
        hit = selectedIdx>=first & selectedIdx<=final;
        if any(hit)
            blockSums(end+1,1)=sum(amplitude(hit).*exp(1i*phase(selectedIdx(hit)))); %#ok<AGROW>
        end
        first=final+1;
    end
end
nBlocks=numel(blockSums);
if nBlocks<4
    error('PACT:InsufficientBlocks', ...
        'Only %d HAS-containing blocks; at least four are required for exploratory surrogate inference.',nBlocks);
end

miNull=zeros(1,nSurro);
nHAS=numel(selectedIdx);
% Chunk surrogates to cap peak memory, and use the caller's RNG state.
chunkSize=256;
for first=1:chunkSize:nSurro
    last=min(first+chunkSize-1,nSurro);
    count=last-first+1;
    phaseOffsets=2*pi*rand(nBlocks,count);
    miNull(first:last)=abs(sum(bsxfun(@times,blockSums,exp(1i*phaseOffsets)),1))/nHAS;
end
end
