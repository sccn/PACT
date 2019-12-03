% pac_PhaseDepVis() 
%
% Usage:
%   >> EEG = pac_PhaseDepVis(EEG);
%
% Inputs:
%  EEG     : EEGLAB structure
%
% Dependency:
%  This function calls myPolar.m which is contained in this plugin.
%
% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD 2012-
% 
% 10/31/2012 ver 2.0 by Makoto. Phase-sorted amp updated.
% 10/22/2012 ver 1.0 by Makoto.

% Copyright (C) 2012 Makoto Miyakoshi, JSPS/SCCN,INC,UCSD
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

function pac_PhaseDepVis(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% select histogram or polar plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    userInput = inputgui('title', 'pac_PhaseDepVis()', 'geom',...
        {{1 3 [0 0] [1 1]} ...
         {1 3 [0 2] [1 1]}},...
        'uilist',...
        {{'style' 'popupmenu' 'string' 'Probability density (histogram)|Probability density (polar plot)|Phase-sorted amp' 'tag' 'phaseplot' 'value' 1} ...
         {'style' 'text'      'string' 'Significant values of modulation index (instantaneous) are shown.'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate and visualize phase dependency %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'Name', 'pac_PhaseDepVis()', 'NumberTitle','off')
set(gcf,'Color', [0.93 0.96 1]);
counter    = 0;
figNumber  = 1;
maxFigNum = ceil(EEG.nbchan/30);
annotation(gcf,'textbox', [0 1 1 0],...
    'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase (page ' num2str(figNumber) ' of ' num2str(maxFigNum) ')']},...
    'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');


if     userInput{1,1} == 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot histogram of probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for n = 1:EEG.nbchan
        if counter == 30
            % this figure is now full - wrap it up and open up a new one
            axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');
            set(gcf, 'Name', 'pac_PhaseDepVis()', 'NumberTitle','off');
            figNumber  = figNumber+1;
            figure
            set(gcf,'Color', [0.93 0.96 1]);
            annotation(gcf,'textbox', [0 1 1 0],...
                'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase (page ' num2str(figNumber) ' of ' num2str(maxFigNum) ')']},...
                'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');

            counter = 0;
        end
        counter = counter+1;

        % open subplot
        subplot(6,5,counter)

        % visualize hist plot
        phasetmp     = EEG.pac.pacAngle{n,1};
        phasetmpPosi = phasetmp(phasetmp>0);
        phasetmpNega = phasetmp(phasetmp<0);
        phasetmpNega = 2*pi+phasetmpNega;
        phasetmp     = [phasetmpPosi phasetmpNega];
        hist(phasetmp, 36);
        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', [16/255 67/255 100/255]); % Marine Blue
        hold on;
        
        % overplot a cosine wave
        t = 0:0.01:2*pi;
        plot(t, (cos(t)+1)/4*max(EEG.pac.histBinMax), 'r', 'LineWidth',2)
        set(gca, 'XTickLabel',{'0','pi','2*pi'},'XTick',[0 3.1415 6.2832], 'CLim',[1 2])
        xlim([0 2*pi+0.01])
        ylim([0 max(EEG.pac.histBinMax)*1.01])
        xlabel('phase')
        ylabel(['probability*' num2str(length(phasetmp))])
        title(['Ch' num2str(n)])
    end
    clear N X phasetmp* t


    
elseif userInput{1,1} == 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot polarplot of probability %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for n = 1:EEG.nbchan
        if counter == 30
            % this figure is now full - wrap it up and open up a new one
            axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');
            set(gcf, 'Name', 'pac_PhaseDepVis()', 'numbertitle','off');
            figNumber  = figNumber+1;
            figure
            set(gcf,'Color', [0.93 0.96 1]);
            annotation(gcf,'textbox', [0 1 1 0],...
                'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase (page ' num2str(figNumber) ' of ' num2str(maxFigNum) ')']},...
                'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');

            counter = 0;
        end
        counter = counter+1;
        
        % open subplot
        subplot(6,5,counter)

        % visualize polar plot
        [N, X] = hist(EEG.pac.pacAngle{n,1}, 36);
        THETA = cat(2, X, X(1));
        RHO = cat(2, N, N(1));
        
        % polar(THETA, RHO);
        myPolar(THETA, RHO, '-', max(EEG.pac.histBinMax)*1.01, 3);
        title(['Ch' num2str(n)])
    end
    clear N X THETA RHO


    
elseif userInput{1,1} == 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot histogram of amplitude %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for n = 1:EEG.nbchan
        if counter == 30
            % this figure is now full - wrap it up and open up a new one
            axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');
            set(gcf, 'Name', 'pac_PhaseDepVis()', 'numbertitle','off');
            figNumber  = figNumber+1;
            figure
            set(gcf,'Color', [0.93 0.96 1]);
            annotation(gcf,'textbox', [0 1 1 0],...
                'String', {[num2str(EEG.pac.hfoAmp(1)) '-' num2str(EEG.pac.hfoAmp(2)) ' Hz amp coupling to  ' num2str(EEG.pac.lfoPhase(1)) '-' num2str(EEG.pac.lfoPhase(2)) ' Hz phase (page ' num2str(figNumber) ' of ' num2str(maxFigNum) ')']},...
                'HorizontalAlignment','center', 'FontSize',12, 'FitBoxToText','off', 'LineStyle','none');
            counter = 0;
        end
        counter = counter+1;

        % open subplot
        subplot(6,5,counter)
        
        % visualize hist plot with pacAmp
        binStep = 2*pi/36; % 360 / 36 = bin is every 10 degrees
        [angleinfo sortindex] = sort(EEG.pac.pacAngle{n,1});
        
        for binNumber = 1:36
            [dummy idy] = find(-pi+(binNumber-1)*binStep < angleinfo & angleinfo < -pi+binNumber*binStep); 
            tmpIndex = sortindex(idy);
            tmpAmp   = EEG.pac.pacAmp{n,1}(1,tmpIndex);
            binAmp(binNumber, 1) = mean(tmpAmp);
            binAmp(binNumber, 2) = std(tmpAmp);
        end
        
        % sort from 0 to 2pi
        binAmp = binAmp([19:36 1:18],:);
        
        bar(binAmp(:,1), 'FaceColor', [0 86/255 69/255], 'BarWidth', 1); % Billiard 
        hold on
        errorbar(1:36, binAmp(:,1), zeros(1,36), binAmp(:,2), 'LineStyle','none','Color',[0 0 0])
        xlim([0.5 36.5])
        ylim([0 max(binAmp(:,1)+binAmp(:,2))*1.05])
        xlabel('phase')
        ylabel('amp (uV)')
        hold on
        
        % overplot a cosine wave
        t = 0.5:0.1:36;
        plot(t, (cos(t/5.7)+1)/2*max(binAmp(:,1)+binAmp(:,2))/2, 'r', 'LineWidth',2)
        set(gca, 'XTickLabel',{'0','pi','2*pi'},'XTick',[0 18 36], 'CLim',[1 2])
        title(['Ch' num2str(n)])
    end
end
axcopy(gcf, 'if ~isempty(get(gca, ''''userdata'''')), eval(get(gca, ''''userdata'''')); end');
set(gcf, 'Name', 'pac_PhaseDepVis()', 'numbertitle','off');