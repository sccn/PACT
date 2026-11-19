function EEG = pac_pop_statsSetUp(EEG)
%PAC_POP_STATSSETUP Configure corrected exploratory MI significance plots.
% Historical phase/phase-binned-amplitude tests are disabled because
% temporal autocorrelation and HAS selection invalidate their p-values.
%
% Copyright (C) 2012-2026 Makoto Miyakoshi. GPL-2.0-or-later.

if ~isfield(EEG,'pac') || ~isfield(EEG.pac,'uncorrected')
    error('Compute PACT before configuring statistics.');
end
methods={'uncorrected','Bonferroni','BonfHolm','FDR'};
if isfield(EEG.pac,'multiCompType') && ismember(EEG.pac.multiCompType,1:4)
    initialMethod=EEG.pac.multiCompType;
else
    initialMethod=3;
end
if isfield(EEG.pac,'confIntType') && ismember(EEG.pac.confIntType,[1 2])
    initialCI=EEG.pac.confIntType;
else
    initialCI=1;
end
choice=inputgui('title','PACT: exploratory MI statistics',...
    'geom',{{2 4 [0 0] [1 1]} {2 4 [1 0] [1 1]} ...
            {2 4 [0 1] [1 1]} {2 4 [1 1] [1 1]} ...
            {2 4 [0 2] [1 1]} {2 4 [1 2] [1 1]}},...
    'uilist',{{'style' 'text' 'string' 'Significance threshold (0 < p < 1)'} ...
               {'style' 'edit' 'string' num2str(EEG.pac.alpha)} ...
               {'style' 'text' 'string' 'Channel-wise multiple comparisons'} ...
               {'style' 'popupmenu' 'string' 'Uncorrected|Bonferroni|Holm|Benjamini-Hochberg FDR' 'value' initialMethod} ...
               {'style' 'text' 'string' 'Surrogate MI null percentile'} ...
               {'style' 'popupmenu' 'string' '95%|99%' 'value' initialCI}});
if isempty(choice)
    return
end
alpha=pac_parseNumeric(choice{1});
if ~isscalar(alpha) || alpha<=0 || alpha>=1
    error('Significance threshold must be strictly between 0 and 1.');
end
selected=choice{2};
if ~isscalar(selected) || ~ismember(selected,1:4)
    error('Invalid multiple-comparison selection.');
end
EEG.pac.alpha=alpha;
EEG.pac.multiCompType=selected;
EEG.pac.confIntType=choice{3};
EEG.pac.currentMIpval=EEG.pac.(methods{selected}).MIpval;
EEG.pac.currentPhaseTest=nan(EEG.nbchan,1);
EEG.pac.currentAmpTest=nan(EEG.nbchan,1);
EEG.pac.currentWtsnWillPval=nan(EEG.nbchan,1);
EEG.pac.currentPhaseTestName='Not tested';
EEG.pac.currentAmpTestName='Not tested';
EEG.pac.currentMcompName=methods{selected};
if choice{3}==1
    EEG.pac.currentModIndConfInt='95%';
else
    EEG.pac.currentModIndConfInt='99%';
end
fprintf(['PACT: MI significance uses a temporally structured surrogate null. ',...
         'Circular/phase-sorted amplitude distribution p-values are disabled.\n']);
end
