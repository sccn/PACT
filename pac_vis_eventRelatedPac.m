function pac_vis_eventRelatedPac(EEG)
%PAC_VIS_EVENTRELATEDPAC Plot event-locked HFO amplitude and circular phase.
% Mean phases are circular means; phase wrap makes additive +/- SD shading
% misleading, so only the mean phase trace is plotted.
%
% Copyright (C) 2024-2026 Makoto Miyakoshi, Henrico Stam.
% GPL-2.0-or-later.

if ~isfield(EEG,'pac') || ~isfield(EEG.pac,'windowMeanAmp')
    warndlg('Event-related PACT data have not been calculated.');
    return
end
answer=inputgui('title','Event-related PACT','geometry',{[1 1]},...
    'uilist',{{'style','text','string','Channel number'},...
              {'style','edit','string','1'}});
if isempty(answer)
    return
end
ch=str2double(answer{1});
if ~isscalar(ch) || ~isfinite(ch) || ch~=round(ch) || ch<1 || ch>EEG.nbchan
    error('Channel number must be an integer between 1 and EEG.nbchan.');
end
meanAmp=EEG.pac.windowMeanAmp(ch,:);
sdAmp=EEG.pac.windowStdAmp(ch,:);
meanPhase=EEG.pac.windowMeanPhase(ch,:);
n=numel(meanAmp);
timeMs=((1:n)-ceil(n/2))*1000/EEG.srate;

figure('NumberTitle','off','Name','PACT event-related exploratory summary',...
       'Color',[1 1 1]);
yyaxis left
hPhase=plot(timeMs,meanPhase,'b','LineWidth',1.5);
ylabel('Circular mean LFO phase (radians)');
ylim([-pi pi]);
yticks([-pi 0 pi]);
yticklabels({'-\pi','0','\pi'});
yyaxis right
hold on
fill([timeMs fliplr(timeMs)],...
     [meanAmp+sdAmp fliplr(meanAmp-sdAmp)],...
     [1 0 0],'EdgeColor','none','FaceAlpha',0.12);
hAmp=plot(timeMs,meanAmp,'r','LineWidth',1.5);
ylabel('HFO envelope (\muV)');
xlabel('Time from event (ms)');
xline(0,':k');
legend([hPhase hAmp],{'Circular mean phase','Mean HFO amplitude +/- SD'},...
       'Location','best');
title(sprintf('Channel %d, %d events, pooled MI %.3g',ch,...
    size(EEG.pac.winEdgesAfterRejection,1),EEG.pac.mi(ch)));
end
