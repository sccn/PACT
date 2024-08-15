% eegplugin_pac(): A plugin for EEGLAB to provide a solution for computing
%                  phase-amplitude coupling.

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD
% History
% 08/12/2024 Makoto. Event-related PAC visualization supported. 
% 01/13/2021 Makoto. Mobilab express removed because unable to update the interactive graphics.
% 07/24/2019 Makoto. scanHfoByPhase added.
% 03/28/2014 ver 1.5 by Makoto. pac_setPath added.
% 03/27/2013 ver 1.4 by Makoto. pac_setMobilabExpress added.
% 01/31/2013 ver 1.3 by Makoto. GUI menu name is officially 'PACT'
% 01/16/2013 ver 1.2 by Makoto. New menu added.
% 12/26/2012 ver 1.1 by Makoto. Minor changes made.
% 10/22/2012 ver 1.0 by Makoto. Created.

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

function vers = eegplugin_pac(fig, try_strings, catch_strings)

% sanity check
vers = 'PACT0.60';
if nargin < 3
    error('eegplugin_pac requires 3 arguments');
end

% run PACT initial path set up
ok = pac_setPath;
if ~ok
    fprintf('PACT initialization failed!\n');
    pause(2);
    return
end

% create a highLevelMenu
highLevelManu = findobj(fig, 'tag', 'tools');
submenu       = uimenu(highLevelManu, 'label', 'PACT','separator','on');

% add submenu
% uimenu( submenu, 'label', 'Invert polarity',                 'callback', 'EEG = pac_pop_invertPolarity(EEG);  [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG); eeglab redraw');
% uimenu( submenu, 'label', 'Convert event names to strings',  'callback', 'EEG = pac_pop_convertNumEventType2Str(EEG);  [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG,EEG); eeglab redraw');
% uimenu( submenu, 'label', 'Handpick HFO (Mobilab)',          'callback', 'pac_setMobilabExpress; pac_plotEdit(EEG);');
% uimenu( submenu, 'label', 'Count and save HFO markings',     'callback', 'pac_pop_countHfo(EEG);');
uimenu( submenu, 'label', 'Compute PAC',                     'callback', '[ALLEEG EEG] = pac_pop_main(ALLEEG, EEG, CURRENTSET);','separator','on');
uimenu( submenu, 'label', 'Plot HFO-marked raw data',        'callback', 'pac_vis_rawDataPlot(EEG);');
uimenu( submenu, 'label', 'Copy event markers',              'callback', '[ALLEEG,EEG] = pop_copyEventMarker(ALLEEG,EEG);');
uimenu( submenu, 'label', 'Set up statistics',               'callback', 'EEG = pac_pop_statsSetUp(EEG);','separator','on');
uimenu( submenu, 'label', 'Plot Modulation Index',           'callback', 'pac_vis_chanMi(EEG);','separator','on');
uimenu( submenu, 'label', 'Plot angular hist (bar)',         'callback', 'pac_vis_angHistBar(EEG);');
uimenu( submenu, 'label', 'Plot angular hist (polar)',       'callback', 'pac_vis_angHistPolar(EEG);');
uimenu( submenu, 'label', 'Plot phase-sorted amp',           'callback', 'pac_vis_phaseSortAmp(EEG);');
uimenu( submenu, 'label', 'Plot event-related PAC',          'callback', 'pac_vis_eventRelatedPac(EEG);');
uimenu( submenu, 'label', 'Scan HAS x Phase_Freq',           'callback', 'EEG = pac_pop_scanLfoPhaseFreq(EEG);','separator','on');
uimenu( submenu, 'label', 'Scan HFO_Freq x Phase_Freq',      'callback', 'EEG = pac_pop_scanHfoByPhase(EEG);');
%uimenu( submenu, 'label', 'Scan LFO freqs (EZ, with phase)', 'callback', 'pac_pop_compareLfoFreqs(EEG);');
  