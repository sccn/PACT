% pac_man() - Computes phase-amplitude coupling measure. Requires inputs
%             from pac_pop_main. 
% Usage:
%   >> EEG = pac_man(EEG, varargin);
%   >> EEG = pac_man(EEG, 'lfoPhase', [1,4], 'hfoAmp', [80,300],...
%             'hfoTopRatio', 10, 'hfoPool', 1, 'whichMarker', 1,...
%             'windowLength', [], 'alpha', 0.05, 'numSurro', 500,...
%             'numPhaseBin', 30);
%
% Input arguments (ALL necessary):
%     lfoPhase     : Band-pass frequencies for low-frequency phase.
%     hfoAmp       : Band-pass frequencies for high-frequency oscillation
%                    (HFO) amplitude.
%     hfoTopRatio  : Highest amplitude rate in HFO (unit: percent)
%     hfoPool      : 1, each channel; 2, all channels; 3, use user-selected
%                    event markers. Use 1 for the optimal automatic processes.
%     whichMarker  : When user-selected event marker is specified for 'hfoPool',
%                    select the type of event to use for event markers.
%                    Otherwise, can be left empty. The index of the event
%                    marker used can be identified by running 
%                    eventTypes = unique({EEG.event.type}).
%     windowLength : When user-selected event marker is specified for 'hfoPool',
%                    select the window length. Ex. 500 means a window of
%                    +/- 500 ms from the center of the event marker will be
%                    used for a window. Overlapped winodws are not counted
%                    only once.         
%     alpha        : statistical significance level. Eg. 0.05 -> 5% level.
%     numSurro     : Number of iteration for surrogation test. This test
%                    method is explained in Canolty et al. (2006). For a
%                    preliminary test, use 200-500. For publication, use 2000-.
%     numPhaseBin  : Number of equidistant phase bins. Increasing this
%                    number results in increase of statistical power. Use
%                    20-30 for start.
%
% Outputs: All outputs are located in EEG.pac.
%     analyticEEG      : Band-pass filterd low-frequency phase coupled with
%                        band-pass filtered high-frequency oscillation
%                        amplitude. Same size as EEG.data.
%     hfoIndex         : Data indices of selected highest amplitudes.
%     angleMean        : Mean angles of data selected with hfoIndex.
%     vectLength       : Resultant vector length of data selected with hfoIndex.
%     vectLengthVar    : Variance of the resultant vector length.
%     histBinMax       : Maximum values in histograms for drawing graphs.
%     phaseProbability : Number of data count in each phase bins.
%     phaseSrotedAmp   : Mean amplitdes of each phase bins.
%     phaseSortedAmpSe : Standard errors of amplitdes in each phase bins.
% 
% These four structures contains results from different multiple comparison corrections. 
%     uncorrected      : Channel-wise multiple comparison uncorrected.
%     Bonferroni       : Channel-wise multiple comparison with Bonferroni.
%     BonfHolm         : Channel-wise multiple comparison with Bonferroni-Holm.
%     FDR              : Channel-wise multiple comparison with false discovery rate.
%
% Each of the structure contains p-values from various tests:
%     MIpval            : Results from phase permutation on Modulation Index.
%     phaseRayleighPval : Results from Rayleigh test on phase distribution.   
%     phaseOmniTestPval : Results from Omnibus test on phase distribution.
%     phaseRaoSpacePval : Results from Rao's spacing test on phase distribution.
%                         This is currently not used for weired results.
%     phaseChi2GofPval  : Results from Chi-square goodness of fit test on phase distribution.
%     phaseKstestPval   : Results from Kolmogorov-Smirnov test on phase distribution.
%     ampChi2GofPval    : Results from Chi-square goodness of fit test on phase-bin-mean amplitude across phase bins.
%     ampKstestPval     : Results from Kolmogorov-Smirnov test on phase-bin-mean amplitude across phase bins.
%
% In the MIpval, the observed Modulation Index is tested against the null
% hypothesis that there is no phase-amplitude coupling in the observed data.
% In the phaseXXXPval, the phase distribution of the HFO-indexed data are tested
% against the null hypothesis that the observed phase distributions are uniform.
% In the ampXXXPval, the phase-bin-mean amplitudes are tested against the
% null hypothesis that the observed phase-bin-mean amplitudes are uniformly
% distributed.
%
% Note that all of these statistics except for MIpval are critically affected
% by numPhaseBin. If you set numPhaseBin to 10000 (suppose your data is longer
% than that), p-value will be ridiculously small (such as 1.616199e-35).

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD 2012.
% BGM RDJ Pacman Power-Pill Mix.
%
% History:
% 07/24/2019 Makoto. 'binValue' size mismatch issue fixed again. Thanks Brian Kavanaugh!
% 05/10/2019 Makoto. 'binValue' size mismatch issue fixed. Thanks Brian Kavanaugh!
% 02/01/2017 ver 5.2 by Makoto. Abnormally low amplitude channel error added.
% 06/03/2014 ver 5.1 by Makoto. phasebin+1 changed. 
% 06/21/2013 ver 5.0 by Makoto. Watson-Williams test added.
% 04/16/2013 ver 4.1 by Makoto. Confidence Interval option added.
% 12/26/2012 ver 4.0 by Makoto. Optimized for Matrix calculation. Statistics upgraded. Help written.
% 11/02/2012 ver 3.0 by Makoto. Statistics and AA added.
% 10/31/2012 ver 2.1 by Makoto. fixed EEG.pac.maskESC_N = pacAllChan4 < critical4N;
% 10/22/2012 ver 2.0 by Makoto. Made into a part of plugin.
% 09/18/2012 ver 1.1 by Makoto. More accurate window length. EEG.pac.
% 09/14/2012 ver 1.0 by Makoto. Prototype created.

% Copyright (C) 2012 Makoto Miyakoshi, JSPS/SCCN,INC,UCSD;
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function EEG = pac_man(EEG, varargin)

%%%%%%%%%%%%%%%%%%%
%%% input check %%%
%%%%%%%%%%%%%%%%%%%
if nargin < 18
    help pac_man
    return
end

result = finputcheck(varargin, {...
    'lfoPhase'     'real'    [0 Inf] []
    'hfoAmp'       'real'    [0 Inf] []
    'hfoTopRatio'  'real'    [0 100] []
    'hfoPool'      'integer' [1 3]   []
    'whichMarker'  'integer' [1 Inf] []
    'windowLength' 'integer' [1 Inf] []
    'alpha'        'real'    [0 100] []
    'numSurro'     'integer' [1 Inf] []
    'numPhaseBin'  'integer' [1 Inf] []});
if ischar(result), error(result); end;
clear result

% Check dead channels
allChanStd = std(EEG.data, [], 2);
allChanStdRatio = allChanStd/median(allChanStd);
if any(allChanStdRatio<0.1)
    badChanIdx = find(allChanStdRatio<0.1);
    for chIdx = 1:length(badChanIdx)
        disp(sprintf('Channel %.0f have abnormally low SD. Remove this channel from analysis.', badChanIdx(chIdx)))
    end
    error('Use ''Edit''->''Select data'' to remove the channels above.')
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove existing pac %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(EEG, 'pac')
    EEG = rmfield(EEG, 'pac');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial parameter setting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfoPhase     = vararginReader('lfoPhase', varargin);
hfoAmp       = vararginReader('hfoAmp', varargin);
hfoTopRatio  = vararginReader('hfoTopRatio', varargin);
hfoPool      = vararginReader('hfoPool', varargin);
whichMarker  = vararginReader('whichMarker', varargin);
windowLength = vararginReader('windowLength', varargin);
alpha        = vararginReader('alpha', varargin);
numSurro     = vararginReader('numSurro', varargin);
numPhaseBin  = vararginReader('numPhaseBin', varargin);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% store parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%
EEG.pac.lfoPhase     = lfoPhase;
EEG.pac.hfoAmp       = hfoAmp;
EEG.pac.hfoTopRatio  = hfoTopRatio;
EEG.pac.hfoPool      = hfoPool;
EEG.pac.whichMarker  = whichMarker;
EEG.pac.windowLength = windowLength;
EEG.pac.alpha        = alpha;
EEG.pac.numSurro     = numSurro;
EEG.pac.numPhaseBin  = numPhaseBin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reserve memory for storing results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEG.pac.analyticEEG  = zeros(size(EEG.data));

%%%%%%%%%%%%%%%%%%
%%% DC removal %%%
%%%%%%%%%%%%%%%%%%    
EEG.data = EEG.data - repmat(mean(EEG.data, 2), [1 length(EEG.data)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute analytic EEG %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% low frequency phase
EEG_lowFreqPhase  = pop_eegfiltnew(EEG, lfoPhase(1), lfoPhase(end));
lowFreqPhase_data = EEG_lowFreqPhase.data;
analyPhase        = angle(hilbert(lowFreqPhase_data'))';

% high frequency amp
EEG_highFreqAmp   = pop_eegfiltnew(EEG, hfoAmp(1),   hfoAmp(end));
highFreqAmp_data  = EEG_highFreqAmp.data;
analyAmp          = abs(hilbert(highFreqAmp_data'))';

% low freq phase + high freq amp
tmpZ = analyAmp.*exp(1i*analyPhase);
EEG.pac.analyticEEG = tmpZ;

clear EEG_* tmp* z_*

switch hfoPool
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% channel-wise pooling %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    case 1
        analyAmpSort   = sort(analyAmp, 2, 'descend');
        critical       = floor(length(analyAmpSort)*hfoTopRatio/100);
        criticalValues = analyAmpSort(:, critical);
        logicalHasMap = analyAmp >= repmat(criticalValues, [1 size(analyAmp,2)]);
        for chIdx = 1:EEG.nbchan
            EEG.pac.hfoIndex{chIdx,1} = find(logicalHasMap(chIdx,:));
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% whole-data pooling %%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 
    case 2
        tmp = analyAmp(:);
        tmpSort = sort(tmp, 'descend');
        critical = tmpSort(floor(length(tmpSort)*EEG.pac.hfoTopRatio/100));
        globalHfoIndex = tmp >= critical;

        % obtain globalHFOmask
        globalHFOmask = false(length(tmpSort),1);
        globalHFOmask(globalHfoIndex) = 1;
        globalHFOmask = reshape(globalHFOmask, size(EEG.data));
        for chIdx = 1:EEG.nbchan
            EEG.pac.hfoIndex{chIdx,1} = find(globalHFOmask(chIdx,:));
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% handpicked window centers +/- halfWinLen  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3 % handpicked
        % create a list that has channel and latency of the selected event 
        eventTypes = unique({EEG.event.type});
        tmpMarker = eventTypes{whichMarker-1};
        hitIndex = strcmp({EEG.event.type}, tmpMarker);
        hitIndex = find(hitIndex);
        tmpChanList    = cell2mat({EEG.event(1,hitIndex).channel})';
        tmpLatencyList = cell2mat({EEG.event(1,hitIndex).latency})';
        hfoList        = cat(2, tmpChanList, tmpLatencyList);
        
        % create HFO index
        EEG.pac.hfoIndex{EEG.nbchan,1} = []; % prepare cell array
        for chIdx = 1:EEG.nbchan
            tmpHfoListIndex = hfoList(:,1)==chIdx;
            if any(tmpHfoListIndex)
                tmpHfoLatencyList = hfoList(tmpHfoListIndex, 2);
                tmpHfoListWin     = zeros(1,length(tmpHfoLatencyList)*windowLength*2);
                for m = 1:length(tmpHfoLatencyList)
                    tmpIndex = (m-1)*windowLength*2+1:m*windowLength*2;
                    tmpHfoListWin(1,tmpIndex) = tmpHfoLatencyList(m)-windowLength:tmpHfoLatencyList(m)+windowLength-1;
                end
                tmpHfoListWin(tmpHfoListWin<1)       =[]; % exclude negative latency
                tmpHfoListWin(tmpHfoListWin>EEG.pnts)=[]; % exclude out of range latency
                tmpHfoListWin = unique(tmpHfoListWin);    % exclude overlap
                EEG.pac.hfoIndex{chIdx,1} = tmpHfoListWin;
            end
        end
        
        % create non-empty channel list
        nonEmptyChannelList = false(1,length(EEG.pac.hfoIndex));
        for chIdx = 1:length(EEG.pac.hfoIndex)
             nonEmptyChannelList(chIdx,1) = ~isempty(EEG.pac.hfoIndex{chIdx,1});
        end
end
clear logical* tmp* critical* global* analy* x_* whichMarker windowLength eventTypes hfoList hitIndex

%%%%%%%%%%%%%%%%%%%%%%%
%%% reserve results %%%
%%%%%%%%%%%%%%%%%%%%%%%
EEG.pac.angleMean                     = zeros(EEG.nbchan,1);
EEG.pac.angleStd                      = zeros(EEG.nbchan,1);
EEG.pac.vectLength                    = zeros(EEG.nbchan,1);
EEG.pac.vectLengthVar                 = zeros(EEG.nbchan,1);
EEG.pac.histBinMax                    = zeros(EEG.nbchan,1);
EEG.pac.mi                            = zeros(EEG.nbchan,1);
EEG.pac.moduInd95CI                   = zeros(EEG.nbchan,1);
EEG.pac.moduInd99CI                   = zeros(EEG.nbchan,1);
EEG.pac.phaseBinValues                = zeros(EEG.nbchan, numPhaseBin+1); % 06/03/2014

EEG.pac.uncorrected.MIpval            = ones(EEG.nbchan,1);
EEG.pac.uncorrected.MIpval            = ones(EEG.nbchan,1);
EEG.pac.uncorrected.phaseRayleighPval = ones(EEG.nbchan,1);
EEG.pac.uncorrected.phaseOmniTestPval = ones(EEG.nbchan,1);
EEG.pac.uncorrected.phaseRaoSpacePval = ones(EEG.nbchan,1);
EEG.pac.uncorrected.phaseWtsnWillPval = ones(EEG.nbchan,1);
EEG.pac.uncorrected.ampChi2GofPval    = ones(EEG.nbchan,1);
EEG.pac.uncorrected.ampKstestPval     = ones(EEG.nbchan,1);

EEG.pac.Bonferroni.MIpval             = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.MIpval             = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.phaseRayleighPval  = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.phaseOmniTestPval  = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.phaseRaoSpacePval  = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.phaseWtsnWillPval  = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.ampChi2GofPval     = ones(EEG.nbchan,1);
EEG.pac.Bonferroni.ampKstestPval      = ones(EEG.nbchan,1);

EEG.pac.BonfHolm.MIpval               = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.MIpval               = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.phaseRayleighPval    = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.phaseOmniTestPval    = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.phaseRaoSpacePval    = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.phaseWtsnWillPval    = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.ampChi2GofPval       = ones(EEG.nbchan,1);
EEG.pac.BonfHolm.ampKstestPval        = ones(EEG.nbchan,1);

EEG.pac.FDR.MIpval                    = ones(EEG.nbchan,1);
EEG.pac.FDR.MIpval                    = ones(EEG.nbchan,1);
EEG.pac.FDR.phaseRayleighPval         = ones(EEG.nbchan,1);
EEG.pac.FDR.phaseOmniTestPval         = ones(EEG.nbchan,1);
EEG.pac.FDR.phaseRaoSpacePval         = ones(EEG.nbchan,1);
EEG.pac.FDR.phaseWtsnWillPval         = ones(EEG.nbchan,1);
EEG.pac.FDR.ampChi2GofPval            = ones(EEG.nbchan,1);
EEG.pac.FDR.ampKstestPval             = ones(EEG.nbchan,1);

EEG.pac.phaseProbability              = cell(EEG.nbchan,1);
EEG.pac.phaseSortedAmp                = cell(EEG.nbchan,1);
EEG.pac.phaseSortedAmpSe              = cell(EEG.nbchan,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute PAC and perform statistics  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitBar = waitbar(0,'Computing start');

for chIdx = 1:EEG.nbchan
    
    % update waitbar    
    waitbar(chIdx/EEG.nbchan, waitBar, ['Processing Channel ' num2str(chIdx)]);
    
    hfoIndex = EEG.pac.hfoIndex{chIdx,1};
    if isempty(hfoIndex)
        disp(['Channel ' num2str(chIdx) ' has no HFO marker.'])
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute HFO-indexed data  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find their Angle
        analyAmp   = abs(  EEG.pac.analyticEEG(chIdx,:));
        analyPhase = angle(EEG.pac.analyticEEG(chIdx,:));
        hfoAmp     = analyAmp(hfoIndex);
        hfoPhase   = analyPhase(hfoIndex);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute descriptive statistics %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vector length mean
        EEG.pac.vectLength(chIdx,1)    = circ_r(hfoPhase');
        
        % vector length variance 
        EEG.pac.vectLengthVar(chIdx,1) = circ_var(hfoPhase');    
        
        % angle mean
        tmpAngle    = circ_mean(hfoPhase');
        
        % angle std
        tmpAngleStd = circ_std(hfoPhase');
        
        % converting angle scale from -pi<x<pi to 0<x<2pi
        if tmpAngle < 0;
            tmpAngle = 2*pi+tmpAngle;
        end
        
        % store results
        EEG.pac.angleMean(chIdx,1) = tmpAngle;
        EEG.pac.angleStd(chIdx,1)  = tmpAngleStd;
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute Modulation Index of this HFO %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpMI               = abs(mean(EEG.pac.analyticEEG(chIdx,hfoIndex)));
        EEG.pac.mi(chIdx,1) = tmpMI; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% statistical test for MI %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        minskip = 1;
        maxskip = length(hfoPhase)-1;
        skipIdx = ceil(length(hfoPhase).*rand(numSurro*3,1)); 
        skipIdx(skipIdx>maxskip) = [];
        skipIdx(skipIdx<minskip) = [];
        skipIdx = skipIdx(1:numSurro,1);
        surroPhase  = zeros(length(hfoIndex), numSurro);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% run matrix computation (much faster) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for surroIdx = 1:numSurro
            tmpSurroPhase   = [analyPhase(skipIdx(surroIdx):end) analyPhase(1:skipIdx(surroIdx)-1)]; % randomize phase instead of amp- consider circular shifts
            surroPhase(:,surroIdx) = tmpSurroPhase(hfoIndex)';
        end
        SurroMI = abs(hfoAmp*exp(1i*surroPhase)/length(hfoAmp));
        
        %%%%%%%%%%%%%%%
        %%% p-value %%%
        %%%%%%%%%%%%%%% 
        EEG.pac.uncorrected.MIpval(chIdx,1) = stat_surrogate_pvals(SurroMI, tmpMI, 'upper');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% confidence interval %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        EEG.pac.moduInd95CI(chIdx,1) = prctile(SurroMI, 95);
        EEG.pac.moduInd99CI(chIdx,1) = prctile(SurroMI, 99);

        clear minskip maxskip s skip tmpMI tmpSurroPhase surroPhase sortSurroMI confInt*

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% statistical test for circular distribution of HFO angles %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % 07/23/2019 Makoto. Upon fix request by Brian Kavanaugh.
        sortHfoPhase  = sort(hfoPhase);
        binSize       = length(hfoPhase)/size(EEG.pac.phaseBinValues,2);
        nonIntegerIdx = 1:binSize:length(hfoPhase);
        integerIdx    = round(nonIntegerIdx);
        binValue      = sortHfoPhase(integerIdx);
    
        % % Sample every N points. Note that 'binValue' must contain
        % % sortHfoPhase(1) and sortHfoPhase(end), but the latter one can be
        % % missed by this indexing -- fixed (05/10/2019 Makoto)
        % binValue     = sortHfoPhase(1:binSize:length(hfoPhase));
        % if length(binValue) == size(EEG.pac.phaseBinValues,2)-1
        %     binValue(end+1) = sortHfoPhase(end);
        % end
        
        % % I thought it is this kind of problem, but it is not (07/22/2019)
        % minusPiToPiEdges = linspace(-pi, pi, size(EEG.pac.phaseBinValues,2));
        % [binValue,EDGES] = histcounts(sortHfoPhase, minusPiToPiEdges);

        % binSize      = floor(length(hfoPhase)/numPhaseBin);
        % sortHfoPhase = sort(hfoPhase);
        % binValue     = sortHfoPhase(1:binSize:length(hfoPhase));
        % if length(binValue) == size(EEG.pac.phaseBinValues,2)-1
        %     binValue(end+1) = sortHfoPhase(end);
        % end
        % 
        % % sample every N points
        % binSize      = floor(length(hfoPhase)/numPhaseBin);
        % sortHfoPhase = sort(hfoPhase);
        % binValue     = sortHfoPhase(1:binSize:length(hfoPhase));
        
        % store bin values
        EEG.pac.phaseBinValues(chIdx,:) = binValue'; 
        
        % uniform distribution statistics (3 circular statistics)
        EEG.pac.uncorrected.phaseRayleighPval(chIdx,1) = circ_rtest(binValue);   % Rayleigh test
        EEG.pac.uncorrected.phaseOmniTestPval(chIdx,1) = circ_otest(binValue);   % omnibus or Hodges-Ajne test test
        EEG.pac.uncorrected.phaseRaoSpacePval(chIdx,1) = circ_raotest(binValue); % Rao's spacing test
        
        clear expectedCounts h p st tmpCdf

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% compute phase-sorted amplitude %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        edges = linspace(-pi,pi,numPhaseBin+1);
        [phaseProbability,binIndex] = histc(hfoPhase, edges);
        phaseSortedAmp   = zeros(numPhaseBin,1);
        phaseSortedAmpSe = zeros(numPhaseBin,1);
        for binIdxIdx = 1:numPhaseBin
            phaseSortedAmp(binIdxIdx,1)   = mean(hfoAmp(binIndex==binIdxIdx));
            % phaseSortedAmp(s,1)   = median(hfoAmp(binIndex==s));
            phaseSortedAmpSe(binIdxIdx,1) = std(hfoAmp(binIndex==binIdxIdx))/sqrt(length(find(binIndex==binIdxIdx)));
        end

        % sort from 0 to 2pi
        phaseSortOrder                    = [numPhaseBin/2+1:numPhaseBin 1:numPhaseBin/2];
        EEG.pac.phaseProbability{chIdx,1} = phaseProbability(phaseSortOrder);
        EEG.pac.phaseSortedAmp{chIdx,1}   = phaseSortedAmp(phaseSortOrder);
        EEG.pac.phaseSortedAmpSe{chIdx,1} = phaseSortedAmpSe(phaseSortOrder);

        % statistics (chi-square test)
        edges = linspace(min(phaseSortedAmp),max(phaseSortedAmp),numPhaseBin+1);
        expectedCounts = repmat(numPhaseBin/(length(edges)-1), [1 length(edges)-1]);
        [~,p] = chi2gof(phaseSortedAmp, 'edges', edges, 'expected', expectedCounts);
        EEG.pac.uncorrected.ampChi2GofPval(chIdx,1) = p;

        % statistics (kstest)
        tmpCdf(:,1) = linspace(min(phaseSortedAmp), max(phaseSortedAmp), numPhaseBin);
        tmpCdf(:,2) = linspace(0,1,numPhaseBin);
        [~,p] = kstest(phaseSortedAmp, tmpCdf);
        EEG.pac.uncorrected.ampKstestPval(chIdx,1) = p;     
        
        % keep record of max bin value
        EEG.pac.histBinMax(chIdx,1) = max(phaseProbability);
    end
end

%%%%%%%%%%%%%%%%%%
%%% angle test %%%
%%%%%%%%%%%%%%%%%%
for chIdx = 1:EEG.nbchan
    [EEG.pac.uncorrected.phaseWtsnWillPval(chIdx,1),~] = circ_wwtest(EEG.pac.phaseBinValues(chIdx,:), EEG.pac.angleMean);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% multiple comparisons %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bonferroni
EEG.pac.Bonferroni.MIpval            = EEG.pac.uncorrected.MIpval*EEG.nbchan;
EEG.pac.Bonferroni.phaseRayleighPval = EEG.pac.uncorrected.phaseRayleighPval*EEG.nbchan;
EEG.pac.Bonferroni.phaseOmniTestPval = EEG.pac.uncorrected.phaseOmniTestPval*EEG.nbchan;
EEG.pac.Bonferroni.phaseRaoSpacePval = EEG.pac.uncorrected.phaseRaoSpacePval*EEG.nbchan;
EEG.pac.Bonferroni.phaseWtsnWillPval = EEG.pac.uncorrected.phaseWtsnWillPval*EEG.nbchan;
EEG.pac.Bonferroni.ampChi2GofPval    = EEG.pac.uncorrected.ampChi2GofPval*EEG.nbchan;
EEG.pac.Bonferroni.ampKstestPval     = EEG.pac.uncorrected.ampKstestPval*EEG.nbchan;

% Bonferroni-Holm
EEG.pac.BonfHolm.MIpval            = bonf_holm(EEG.pac.uncorrected.MIpval);
EEG.pac.BonfHolm.phaseRayleighPval = bonf_holm(EEG.pac.uncorrected.phaseRayleighPval);
EEG.pac.BonfHolm.phaseOmniTestPval = bonf_holm(EEG.pac.uncorrected.phaseOmniTestPval);
EEG.pac.BonfHolm.phaseRaoSpacePval = bonf_holm(EEG.pac.uncorrected.phaseRaoSpacePval);
EEG.pac.BonfHolm.phaseWtsnWillPval = bonf_holm(EEG.pac.uncorrected.phaseWtsnWillPval);
EEG.pac.BonfHolm.ampChi2GofPval    = bonf_holm(EEG.pac.uncorrected.ampChi2GofPval);
EEG.pac.BonfHolm.ampKstestPval     = bonf_holm(EEG.pac.uncorrected.ampKstestPval);

% False-discovery Rate
EEG.pac.FDR.MIpval            = fdr(EEG.pac.uncorrected.MIpval);
EEG.pac.FDR.phaseRayleighPval = fdr(EEG.pac.uncorrected.phaseRayleighPval);
EEG.pac.FDR.phaseOmniTestPval = fdr(EEG.pac.uncorrected.phaseOmniTestPval);
EEG.pac.FDR.phaseRaoSpacePval = fdr(EEG.pac.uncorrected.phaseRaoSpacePval);
EEG.pac.FDR.phaseWtsnWillPval = fdr(EEG.pac.uncorrected.phaseWtsnWillPval);
EEG.pac.FDR.ampChi2GofPval    = fdr(EEG.pac.uncorrected.ampChi2GofPval);
EEG.pac.FDR.ampKstestPval     = fdr(EEG.pac.uncorrected.ampKstestPval);

% set default stats options
EEG.pac.phaseTestType = 1;
EEG.pac.multiCompType = 1;

% close waitber
close(waitBar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% display a message on the screen %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
pause(1)

disp('                %%%%%%%%%                ');              
disp('            %%%%         %%%%            ');              
disp('         %%%                 %%%         ');              
disp('       %%%          %%%        %%%       ');              
disp('      %%           %%% %          %      ');         
disp('     %              %%%       %%%%       ');              
disp('    %                    %%%%%           ');              
disp('    %                %%%%                %%%         %%%         %%% ');             
disp('    %              %%%                  %%%%%       %%%%%       %%%%%');             
disp('    %                %%%%                %%%         %%%         %%% ');            
disp('     %                   %%%%%           ');             
disp('      %%                      %%%%       ');             
disp('       %%%                        % ');               
disp('         %%%                  %%%%       ');                
disp('            %%%%         %%%%%           ');              
disp('                %%%%%%%%%                ');
disp(char(10))
disp('         I computed them all!            ');

pause(1)
clc

disp('PAC is computed and stored in EEG.pac')



%%%%%%%%%%%%%%%%%%%%
%%% subfunctions %%%
%%%%%%%%%%%%%%%%%%%%
function output = vararginReader(strings,varargin)
varargin = varargin{1,1};
for n = 1:length(varargin)
    if strcmp(varargin{1,n}, strings)
        output = varargin{1,n+1};
    end
end