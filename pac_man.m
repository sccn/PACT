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
%   * Time shifts LFO phase within uninterrupted valid segments, retaining
%     temporal autocorrelation. The amplitude sample set remains fixed.
%   * Uses (1 + exceedances)/(1 + iterations) for empirical p-values.
%   * Reports adjusted channel-wise Holm and Benjamini-Hochberg p-values.
%   * Historical phase and amp distribution tests are not inferentially
%     valid for temporally dependent, HAS-selected samples and are disabled.
%
% Event-related mode (hfoPool=4) pools HAS samples only from valid event
% windows. Its MI p-value is an exploratory time-shift test, NOT a trial-
% level causal or condition-effect inference.
%
% Copyright (C) 2012-2026 Makoto Miyakoshi. GPL-2.0-or-later.

if nargin < 2 || mod(numel(varargin),2) ~= 0
    error('Specify PACT options as name/value pairs.');
end
if ~isfield(EEG,'trials') || EEG.trials ~= 1 || ndims(EEG.data) ~= 2
    error('PACT currently supports continuous two-dimensional EEG only.');
end
if ~isnumeric(EEG.data) || isempty(EEG.data) || any(~isfinite(EEG.data(:)))
    error('EEG.data must contain finite numeric samples.');
end
if EEG.srate <= 0 || size(EEG.data,2) ~= EEG.pnts || size(EEG.data,1) ~= EEG.nbchan
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
workEEG.data = bsxfun(@minus,double(EEG.data),mean(double(EEG.data),2));
lowEEG = pop_eegfiltnew(workEEG,g.lfoPhase(1),g.lfoPhase(2));
highEEG = pop_eegfiltnew(workEEG,g.hfoAmp(1),g.hfoAmp(2));
phase = angle(hilbert(double(lowEEG.data')))'; 
amp = abs(hilbert(double(highEEG.data')))'; 
clear workEEG lowEEG highEEG

% Original sample-coordinate mask: 1 second near edges and boundaries.
guard = max(1,round(EEG.srate));
validMask = true(1,EEG.pnts);
validMask(1:min(guard,EEG.pnts)) = false;
validMask(max(1,EEG.pnts-guard+1):EEG.pnts) = false;
eventTypes = cell(1,numel(EEG.event));
for k = 1:numel(EEG.event)
    typ = EEG.event(k).type;
    if isnumeric(typ), eventTypes{k} = num2str(typ);
    elseif isstring(typ), eventTypes{k} = char(typ);
    else, eventTypes{k} = typ;
    end
end
boundaryIdx = find(strcmpi(eventTypes,'boundary'));
for k = boundaryIdx
    latency = round(EEG.event(k).latency);
    validMask(max(1,latency-guard):min(EEG.pnts,latency+guard)) = false;
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
        span = centers(k)-halfWin:centers(k)+halfWin;
        if span(1)<1 || span(end)>EEG.pnts || ~all(validMask(span))
            continue
        end
        eventMask(span) = true;
        selectedWindows(end+1,:) = [span(1) span(end)]; %#ok<AGROW>
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

if g.hfoPool == 2
    pooledThreshold=prctile(reshape(amp(:,candidateIdx),[],1),100-g.hfoTopRatio);
end
edges=linspace(-pi,pi,g.numPhaseBin+1);
% Fixed seed for analysis; restore caller RNG even when an error occurs.
priorRng=rng;
restoreRng=onCleanup(@()rng(priorRng)); %#ok<NASGU>
rng(0,'twister');
out.rngSeed=0;

for ch=1:EEG.nbchan
    if g.hfoPool==2
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
    out.angleStd(ch)=sqrt(-2*log(max(out.vectLength(ch),eps)));
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

    [nullMI,~,~]=pac_shiftSurrogate(phase(ch,:),selectedAmp,selectedIdx,...
                                  validMask,g.numSurro,EEG.srate,g.numPhaseBin);
    out.uncorrected.MIpval(ch)=(1+sum(nullMI>=out.mi(ch)))/(g.numSurro+1);
    out.moduInd95CI(ch)=prctile(nullMI,95);
    out.moduInd99CI(ch)=prctile(nullMI,99);
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
out.surrogateMethod='within-segment circular phase shift; HAS amplitudes fixed';
out.inferenceNotice=['Exploratory surrogate MI test. Phase-bin distribution tests ',...
    'are disabled; event-related mode is not a trial-level hypothesis test.'];
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
    out.windowStdPhase=sqrt(-2*log(max(abs(mean(exp(1i*winPhase),3)),eps)));
end
EEG.pac=out;
end
