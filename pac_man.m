function EEG = pac_man(EEG, varargin)
%PAC_MAN Compute highest-amplitude-sampled phase-amplitude coupling (PACT).
% EEG = pac_man(EEG, 'lfoPhase',[1 2], 'hfoAmp',[80 150], ...
%   'hfoTopRatio',2, 'hfoPool',1, 'whichMarker',1, ...
%   'windowLength',[], 'alpha',.05, 'numSurro',2000, 'numPhaseBin',36)
%
% Scientific-safety revision:
%   * Does not modify original EEG.data or EEG.event.
%   * Maintains sample coordinates; excluded samples remain in analyticEEG
%     but are never selected or used to build surrogate distributions.
%   * MI (field mi) keeps its historical definition |mean(a*exp(i*phase))|
%     over HAS samples and is descriptive.
%   * The MI p-value tests the block-debiased MI (field dmi): within each
%     boundary-safe block (>= 5 cycles of the lowest LFO frequency, >= 0.5 s)
%     the block's mean candidate-sample phasor is subtracted before summing
%     a*exp(i*phase) over HAS samples. The null randomly flips the sign of
%     each block's contribution. See pac_blockSurrogate for assumptions.
%     The legacy sample-permutation null was invalid for autocorrelated
%     signals and is not used.
%   * Uses (1 + exceedances)/(1 + iterations) for empirical p-values.
%   * Reports adjusted channel-wise Holm and Benjamini-Hochberg p-values.
%   * Historical phase and amp distribution tests are not inferentially
%     valid for temporally dependent, HAS-selected samples and are disabled.
%
% Event-related mode (hfoPool=4) pools HAS samples only from valid event
% windows. Its MI p-value is an exploratory block sign-flip test, NOT a
% trial-level causal or condition-effect inference. Event-locked waveforms
% can produce small p-values without genuine amplitude modulation.
%
% Copyright (C) 2012-2026 Makoto Miyakoshi. GPL-2.0-or-later.

if nargin < 2 || mod(numel(varargin),2) ~= 0
    error('Specify PACT options as name/value pairs.');
end
if ~isfield(EEG,'trials') || EEG.trials ~= 1 || ndims(EEG.data) ~= 2
    error('PACT currently supports continuous two-dimensional EEG only.');
end
if ~isnumeric(EEG.data) || ~isreal(EEG.data) || isempty(EEG.data) || any(~isfinite(EEG.data(:)))
    error('EEG.data must contain finite numeric samples.');
end
if ~isscalar(EEG.srate) || ~isfinite(EEG.srate) || EEG.srate <= 0 || size(EEG.data,2) ~= EEG.pnts || size(EEG.data,1) ~= EEG.nbchan
    error('Invalid EEGLAB sampling rate or data dimensions.');
end

p = inputParser;
p.FunctionName = 'pac_man';
addParameter(p,'lfoPhase',[],@(x) isnumeric(x) && isreal(x) && numel(x)==2);
addParameter(p,'hfoAmp',[],@(x) isnumeric(x) && isreal(x) && numel(x)==2);
addParameter(p,'hfoTopRatio',[],@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<=100);
addParameter(p,'hfoPool',1,@(x) isnumeric(x) && isscalar(x) && ismember(x,1:4));
addParameter(p,'whichMarker',1,@(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
addParameter(p,'windowLength',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x>0));
addParameter(p,'alpha',.05,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<1);
addParameter(p,'numSurro',2000,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=1 && x==round(x));
addParameter(p,'numPhaseBin',36,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=4 && x==round(x) && mod(x,2)==0);
parse(p,varargin{:});
g = p.Results;
nyquist = EEG.srate/2;
if any(~isfinite(g.lfoPhase)) || g.lfoPhase(1)<=0 || g.lfoPhase(1)>=g.lfoPhase(2) || g.lfoPhase(2)>=nyquist
    error('lfoPhase must be [low high], strictly within (0, Nyquist).');
end
if any(~isfinite(g.hfoAmp)) || g.hfoAmp(1)<=0 || g.hfoAmp(1)>=nyquist || ...
        ~(g.hfoAmp(2)==0 || (g.hfoAmp(2)>g.hfoAmp(1) && g.hfoAmp(2)<nyquist))
    error('hfoAmp must be [low high] or [low 0] for high-pass filtering.');
end
if g.hfoPool == 3
    error('MobiLab handpicked sampling is deprecated and unsupported.');
end
if g.hfoPool == 4 && isempty(g.windowLength)
    error('Event-based sampling needs windowLength in milliseconds.');
end

% Filter a private copy, never the returned EEG data.
workEEG = EEG;
if isfield(workEEG,'pac'), workEEG = rmfield(workEEG,'pac'); end
workEEG.data = bsxfun(@minus,double(EEG.data),mean(double(EEG.data),2));
% Give FIRfilt an unambiguous boundary list even with mixed marker types.
isBoundary=false(1,numel(EEG.event));
for k=1:numel(EEG.event)
    value=EEG.event(k).type;
    isBoundary(k)=isequal(value,-99) || ...
        ((ischar(value) || (isstring(value) && isscalar(value))) && strcmpi(value,'boundary'));
end
workEEG.event=EEG.event(isBoundary);
for k=1:numel(workEEG.event)
    workEEG.event(k).type='boundary';
end
starts=findboundaries(workEEG.event);
cuts=[starts(starts>=1 & starts<=EEG.pnts) EEG.pnts+1];
[lowEEG,~,loFir]=pop_eegfiltnew(workEEG,g.lfoPhase(1),g.lfoPhase(2));
[highEEG,~,hiFir]=pop_eegfiltnew(workEEG,g.hfoAmp(1),g.hfoAmp(2));
guard=max(round(EEG.srate),(max(numel(loFir),numel(hiFir))-1)/2);
phase=zeros(size(workEEG.data));
amp=zeros(size(workEEG.data));
validMask=false(1,EEG.pnts);
for segment=1:numel(cuts)-1
    span=cuts(segment):cuts(segment+1)-1;
    if numel(span)<=2*guard, continue, end
    phase(:,span)=angle(hilbert(double(lowEEG.data(:,span).'))).';
    amp(:,span)=abs(hilbert(double(highEEG.data(:,span).'))).';
    validMask(span(1)+guard:span(end)-guard)=true;
end
clear workEEG lowEEG highEEG
eventTypes = cell(1,numel(EEG.event));
for k = 1:numel(EEG.event)
    typ = EEG.event(k).type;
    if isnumeric(typ), eventTypes{k} = num2str(typ);
    elseif isstring(typ), eventTypes{k} = char(typ);
    else, eventTypes{k} = typ;
    end
end
if ~any(validMask)
    error('No valid samples remain after boundary exclusion.');
end

eventMask = true(1,EEG.pnts);
selectedWindows = zeros(0,2);
if g.hfoPool == 4
    eventMask(:) = false;
    allTypes = unique(eventTypes);
    if g.whichMarker > numel(allTypes)
        error('whichMarker exceeds number of event types.');
    end
    eventIdx = find(strcmp(eventTypes,allTypes{g.whichMarker}));
    halfWin = round(g.windowLength * EEG.srate / 1000);
    if halfWin < 1
        error('Event window is shorter than one sample.');
    end
    centers = round([EEG.event(eventIdx).latency]);
    for k = 1:numel(centers)
        first=centers(k)-halfWin;
        last=centers(k)+halfWin;
        if ~isfinite(centers(k)) || first<1 || last>EEG.pnts
            continue
        end
        span=first:last;
        if all(validMask(span))
            eventMask(span)=true;
            selectedWindows(end+1,:)=[first last]; %#ok<AGROW>
        end
    end
    if isempty(selectedWindows)
        error('No valid event windows remain.');
    end
end
candidateMask = validMask & eventMask;
candidateIdx = find(candidateMask);
if isempty(candidateIdx)
    error('No valid data for PAC.');
end

out = struct;
out.lfoPhase=g.lfoPhase;
out.hfoAmp=g.hfoAmp;
out.hfoTopRatio=g.hfoTopRatio;
out.hfoPool=g.hfoPool;
out.whichMarker=g.whichMarker;
out.windowLength=g.windowLength;
out.alpha=g.alpha;
out.numSurro=g.numSurro;
out.numPhaseBin=g.numPhaseBin;
out.validSampleMask=validMask;
out.selectedSampleMask=candidateMask;
out.analyticEEG=amp.*exp(1i*phase);
out.hfoIndex=cell(EEG.nbchan,1);
out.mi=nan(EEG.nbchan,1);
out.dmi=nan(EEG.nbchan,1);
out.numBlocks=nan(EEG.nbchan,1);
out.blockLength=NaN;
out.vectLength=nan(EEG.nbchan,1);
out.vectLengthVar=nan(EEG.nbchan,1);
out.angleMean=nan(EEG.nbchan,1);
out.angleStd=nan(EEG.nbchan,1);
out.moduInd95CI=nan(EEG.nbchan,1);
out.moduInd99CI=nan(EEG.nbchan,1);
out.phaseBinValues=nan(EEG.nbchan,g.numPhaseBin+1);
out.histBinMax=zeros(EEG.nbchan,1);
out.phaseProbability=cell(EEG.nbchan,1);
out.phaseSortedAmp=cell(EEG.nbchan,1);
out.phaseSortedAmpSe=cell(EEG.nbchan,1);
out.uncorrected.MIpval=nan(EEG.nbchan,1);

if g.hfoPool == 2 && g.hfoTopRatio < 100
    pooledThreshold=prctile(reshape(amp(:,candidateIdx),[],1),100-g.hfoTopRatio);
end
edges=linspace(-pi,pi,g.numPhaseBin+1);
% Fixed seed for analysis; restore caller RNG even when an error occurs.
priorRng=rng;
restoreRng=onCleanup(@()rng(priorRng)); %#ok<NASGU>
rng(0,'twister');
out.rngSeed=0;

for ch=1:EEG.nbchan
    if g.hfoTopRatio==100
        threshold=-Inf;
    elseif g.hfoPool==2
        threshold=pooledThreshold;
    else
        threshold=prctile(amp(ch,candidateIdx),100-g.hfoTopRatio);
    end
    selectedIdx=candidateIdx(amp(ch,candidateIdx)>threshold);
    out.hfoIndex{ch}=selectedIdx;
    if numel(selectedIdx)<max(4,g.numPhaseBin)
        warning('PACT:LowHAS','Channel %d: insufficient HAS samples; statistics set to NaN.',ch);
        continue
    end
    selectedPhase=phase(ch,selectedIdx);
    selectedAmp=amp(ch,selectedIdx);
    phasors=exp(1i*selectedPhase);
    out.mi(ch)=abs(mean(selectedAmp.*phasors));
    out.vectLength(ch)=abs(mean(phasors));
    out.vectLengthVar(ch)=1-out.vectLength(ch);
    out.angleMean(ch)=mod(angle(mean(phasors)),2*pi);
    out.angleStd(ch)=sqrt(2*max(0,1-out.vectLength(ch)));
    sortedPhase=sort(selectedPhase);
    out.phaseBinValues(ch,:)=quantile(sortedPhase,linspace(0,1,g.numPhaseBin+1));
    bins=discretize(selectedPhase,edges);
    counts=accumarray(bins(:),1,[g.numPhaseBin 1],@sum,0);
    means=accumarray(bins(:),selectedAmp(:),[g.numPhaseBin 1],@mean,NaN);
    stds=accumarray(bins(:),selectedAmp(:),[g.numPhaseBin 1],@std,NaN);
    se=stds./sqrt(counts);
    ordering=[g.numPhaseBin/2+1:g.numPhaseBin 1:g.numPhaseBin/2];
    out.phaseProbability{ch}=counts(ordering)';
    out.phaseSortedAmp{ch}=means(ordering);
    out.phaseSortedAmpSe{ch}=se(ordering);
    out.histBinMax(ch)=max(counts);

    try
        [out.dmi(ch),nullStat,info]=pac_blockSurrogate(phase(ch,:),selectedAmp,...
            selectedIdx,validMask,candidateMask,g.numSurro,EEG.srate,g.lfoPhase(1));
        out.numBlocks(ch)=info.nBlocks;
        out.blockLength=info.blockLength;
        out.uncorrected.MIpval(ch)=(1+sum(nullStat>=out.dmi(ch)))/(g.numSurro+1);
        out.moduInd95CI(ch)=prctile(nullStat,95);
        out.moduInd99CI(ch)=prctile(nullStat,99);
    catch ME
        if strcmp(ME.identifier,'PACT:InsufficientBlocks')
            warning('PACT:InsufficientBlocks','Channel %d: %s',ch,ME.message);
        else
            rethrow(ME)
        end
    end
end

% Legacy tests were incorrectly applied to HAS-selected temporal samples.
% Keep field names for plotting compatibility but do not report p-values.
legacy={'phaseRayleighPval','phaseOmniTestPval','phaseRaoSpacePval',...
        'phaseWtsnWillPval','ampChi2GofPval','ampKstestPval'};
for k=1:numel(legacy)
    out.uncorrected.(legacy{k})=nan(EEG.nbchan,1);
end
families=fieldnames(out.uncorrected);
methods={'Bonferroni','BonfHolm','FDR'};
names={'bonferroni','holm','fdr'};
for j=1:numel(methods)
    for k=1:numel(families)
        out.(methods{j}).(families{k})=pac_adjustPValues(...
            out.uncorrected.(families{k}),names{j});
    end
end
out.phaseTestType=1;
out.phaseSortAmpTestType=1;
out.multiCompType=3;
out.confIntType=1;
out.currentMIpval=out.BonfHolm.MIpval;
out.currentPhaseTest=nan(EEG.nbchan,1);
out.currentAmpTest=nan(EEG.nbchan,1);
out.currentWtsnWillPval=nan(EEG.nbchan,1);
out.currentPhaseTestName='descriptive';
out.currentAmpTestName='descriptive';
out.currentMcompName='BonfHolm';
out.currentModIndConfInt='95%';
out.surrogateMethod=['block-debiased MI (dmi); independent random sign flips of ',...
    'nonoverlapping, boundary-safe block contributions'];
out.inferenceNotice=['Exploratory test of HFO amplitude-LFO phase dependence; p-values ',...
    'refer to dmi, not mi. Assumes independent, centrally symmetric block ',...
    'contributions under the null. Phase-bin tests are disabled. Event-locked ',...
    'waveforms and broadband transients can yield small p-values.'];
if g.hfoPool==4
    out.winEdgesAfterRejection=selectedWindows;
    nWin=size(selectedWindows,1);
    winLen=selectedWindows(1,2)-selectedWindows(1,1)+1;
    winAmp=zeros(EEG.nbchan,winLen,nWin);
    winPhase=zeros(EEG.nbchan,winLen,nWin);
    for w=1:nWin
        span=selectedWindows(w,1):selectedWindows(w,2);
        winAmp(:,:,w)=amp(:,span);
        winPhase(:,:,w)=phase(:,span);
    end
    out.windowMeanAmp=mean(winAmp,3);
    out.windowStdAmp=std(winAmp,0,3);
    out.windowMeanPhase=angle(mean(exp(1i*winPhase),3));
    % Angular deviation sqrt(2(1-R)), as CircStat's circ_std and angleStd.
    out.windowStdPhase=sqrt(2*max(0,1-abs(mean(exp(1i*winPhase),3))));
end
EEG.pac=out;
end
