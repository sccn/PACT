% pac_pop_statsSetUp(): Launch GUI to let users select types of statistics
%                       and threshold.
%
% Usage:
%   >>  EEG = pac_pop_statsSetUp(EEG);

% Author: Makoto Miyakoshi, Arnaud Delorme JSPS/SCCN,INC,UCSD
% History
% 04/16/2013 ver 1.2 by Makoto. Confidence Interval option added.
% 12/24/2012 ver 1.1 by Makoto. Minor change added.
% 11/16/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
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

function EEG = pac_pop_statsSetUp(EEG)

try
userInput = inputgui('title', 'pac_pop_statsSetUp()', 'geom', ...
   {{2 5 [0 0] [1 1]}, {2 5 [1 0] [1 1]}...
    {2 5 [0 1] [1 1]}, {2 5 [1 1] [1 1]}...
    {2 5 [0 2] [1 1]}, {2 5 [1 2] [1 1]}...
    {2 5 [0 3] [1 1]}, {2 5 [1 3] [1 1]}...
    {2 5 [0 4] [1 1]}, {2 5 [1 4] [1 1]}},...
    'uilist',...
   {{'style' 'text' 'string' 'Significance threshold [p]'}       {'style' 'edit' 'string' num2str(EEG.pac.alpha)} ...
    {'style' 'text' 'string' 'Test phase distribution with'}     {'style' 'popupmenu' 'string' 'Rayleigh|Omnibus (Hodges-Ajne)|Rao' 'value' EEG.pac.phaseTestType} ...
    {'style' 'text' 'string' 'Test phase-sorted amp with'}       {'style' 'popupmenu' 'string' 'Chi-square goodness of fit|Kolmogorov-Smirnov' 'value' EEG.pac.phaseSortAmpTestType}...
    {'style' 'text' 'string' 'Correct multiple comparison with'} {'style' 'popupmenu' 'string' 'no correction|Bonferroni|Benferroni-Holm|False Discovery Rate' 'value' EEG.pac.multiCompType}...
    {'style' 'text' 'string' 'Confident Interval for MI plot'}   {'style' 'popupmenu' 'string' '95%|99%' 'value' EEG.pac.confIntType}});
catch
userInput = inputgui('title', 'pac_pop_statsSetUp()', 'geom', ...
   {{2 5 [0 0] [1 1]}, {2 5 [1 0] [1 1]}...
    {2 5 [0 1] [1 1]}, {2 5 [1 1] [1 1]}...
    {2 5 [0 2] [1 1]}, {2 5 [1 2] [1 1]}...
    {2 5 [0 3] [1 1]}, {2 5 [1 3] [1 1]}...
    {2 5 [0 4] [1 1]}, {2 5 [1 4] [1 1]}},...
    'uilist',...
   {{'style' 'text' 'string' 'Significance threshold [p]'}       {'style' 'edit' 'string' num2str(EEG.pac.alpha)} ...
    {'style' 'text' 'string' 'Test phase distribution with'}     {'style' 'popupmenu' 'string' 'Rayleigh|Omnibus (Hodges-Ajne)|Rao' 'value' 1} ...
    {'style' 'text' 'string' 'Test phase-sorted amp with'}       {'style' 'popupmenu' 'string' 'Chi-square goodness of fit|Kolmogorov-Smirnov' 'value' 1}...
    {'style' 'text' 'string' 'Correct multiple comparison with'} {'style' 'popupmenu' 'string' 'no correction|Bonferroni|Benferroni-Holm|False Discovery Rate' 'value' 3}...
    {'style' 'text' 'string' 'Confident Interval for MI plot'}   {'style' 'popupmenu' 'string' '95%|99%' 'value' 1}});
end

% store user selection
EEG.pac.alpha                = str2num(userInput{1,1});
EEG.pac.phaseTestType        = userInput{1,2};
EEG.pac.phaseSortAmpTestType = userInput{1,3};
EEG.pac.multiCompType        = userInput{1,4};
EEG.pac.confIntType          = userInput{1,5};

% phase test type
switch userInput{1,2}
    case 1; strPhaseTestType = 'phaseRayleighPval';
    case 2; strPhaseTestType = 'phaseOmniTestPval';
    case 3; strPhaseTestType = 'phaseRaoPval';
end

% phase-sorted amp test type
switch userInput{1,3}
    case 1; strAmpTestType = 'ampChi2GofPval';
    case 2; strAmpTestType = 'ampKstestPval';
end

% multiple comparison correction type
switch userInput{1,4}
    case 1; strMCType = 'uncorrected';
    case 2; strMCType = 'Bonferroni';
    case 3; strMCType = 'BonfHolm';
    case 4; strMCType = 'FDR';
end

% MI confidence interval type
switch userInput{1,5}
    case 1; strMICIType = '95%';
    case 2; strMICIType = '99%';
end

% generate strings
EEG.pac.currentMIpval       = eval(['EEG.pac.' strMCType '.MIpval']);
EEG.pac.currentPhaseTest    = eval(['EEG.pac.' strMCType '.' strPhaseTestType]);
EEG.pac.currentAmpTest      = eval(['EEG.pac.' strMCType '.' strAmpTestType]);
EEG.pac.currentWtsnWillPval = eval(['EEG.pac.' strMCType '.phaseWtsnWillPval']);

EEG.pac.currentPhaseTestName = strPhaseTestType(6:end-4);
EEG.pac.currentAmpTestName   = strAmpTestType(4:end-4);
EEG.pac.currentMcompName     = strMCType;
EEG.pac.currentModIndConfInt = strMICIType;

str = 'Statistics set up successfully.';
disp(str)