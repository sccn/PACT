% pac_plotEdit(): a part of MoBILAB customized for maximally intuitive 
%                EEG event marking. DO NOT USE WITH MOBILAB otherwise paths
%                should conflict.

% History
% 10/08/2015 ver 1.5 by Makoto. Amplitude Scale added.
% 06/03/2014 ver 1.4 by Makoto. Bug fixed in obj.mouseWheelMoveInSec
% 04/28/2014 ver 1.3 by Makoto. 1:obj.streamHandle.numberOfChannels
% 03/29/2014 ver 1.2 by Makoto. Redesigned plotThisTimeStamp() for debugging.
% 03/27/2013 ver 1.1 by Makoto and Alejandro. Mouse scroll. Independent of Mobilab.
% 12/26/2012 ver 1.0 by Alejandro and Makoto. Created.

% Copyright (C) 2013, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
%                     Alejandro Ojeda  SCCN,INC,UCSD
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

%% Creates an object from the class pac_plotEditHandle
%
% pac_plotEditHandle object controls the application pac_plotEdit
%
% Authors: Alejandro Ojeda and Makoto Miyakoshi
% Swartz Center for Computational Neuroscience, Institute for Neural Computation, University
% of California San Diego; Japan Society for the Promotion of Science.

classdef pac_plotEdit < eegplotNGHandle
    properties
        newEventType
        eventMarker
        markerType
        markerSize
        markerEdgeLineWidth
        markerEdgeColor
        cursorType
        cursorSize
        cursorLineWidth
        markerSelectErrorMarginRate
        mouseWheelMoveInSec
     end
    methods
        function obj = pac_plotEdit(EEG,onButtonDownCallback,onMouseMoveCallback)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% shortcut for mouse function %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin < 2, onButtonDownCallback = 'axes_ButtonDownFcn';end
            if nargin < 3, onMouseMoveCallback = 'onMouseMoveCallback';end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% channel status check %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(EEG.chanlocs)
                disp('Channel information is Empty. Generating incremental indices.')
                for n = 1:EEG.nbchan
                    EEG.chanlocs(1,n).labels = num2str(n);
                    EEG.chanlocs(1,n).index  = n;
                end
                
                assignin('base','EEG',EEG);
            end            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% event status check %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(EEG.event)
                disp('Event information is Empty. Generating an initial event.')
                EEG.event.type       = 'start';
                EEG.event.latency    = 1;
                EEG.event.init_index = 1;
                EEG.event.init_time  = 1/EEG.srate;
                EEG.event.urevent = 1;

                EEG.urevent = [];
                EEG.urevent.type       = 'start';
                EEG.urevent.latency    = 1;
                EEG.urevent.init_index = 1;
                EEG.urevent.init_time  = 1/EEG.srate;
                
                assignin('base','EEG',EEG);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% inherit eegplotNGHandle %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            obj@eegplotNGHandle(EEG);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% super small Pointer %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %      1 1 1
            %      1 0 1
            %      1 1 1
            %
            %  the hot spot == 0
            pointerShape           = NaN(16,16);
            pointerShape(1, 1:3)   = 1;
            pointerShape(2, [1 3]) = 1;
            pointerShape(3, 1:3)   = 1;
            pointerHotSpot         = [2,2];
            set(obj.figureHandle,'Pointer','custom', 'PointerShapeCData', pointerShape, 'PointerShapeHotSpot', pointerHotSpot)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% default value setting %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.step                        = obj.windowWidth;
            obj.markerType                  = 'o';
            obj.markerSize                  = 5;
            obj.markerEdgeLineWidth         = 0.5;
            obj.markerEdgeColor             = [0,0,0];
            obj.cursorType                  = '+';
            obj.cursorSize                  = 60;
            obj.cursorLineWidth             = 0.5;
            obj.markerSelectErrorMarginRate = 0.4;
            obj.mouseWheelMoveInSec         = 1;
            
            delete(obj.makeSegmentHandle);
            obj.newEventType = 'HFO'; 
            obj.onscreenDisplay = false;
            set(obj.figureHandle,'WindowButtonMotionFcn', {onMouseMoveCallback,obj},'Name','pac_plotEdit(EEG)');

            if isfield(EEG.event,'channel')
                obj.streamHandle.event = event;
                channel = {EEG.event.channel};
                latency = cell2mat({EEG.event.latency});
                type = {EEG.event.type};
                indHFO = ismember(type,obj.newEventType);
                for it=find(indHFO)
                    if ismac, [~,uuid] = system('uuidgen'); else uuid =  java.util.UUID.randomUUID;end
                    type{it} = [type{it} '_' char(uuid) '_' num2str(channel{it})];
                end
                
                obj.streamHandle.event = obj.streamHandle.event.addEvent(latency(indHFO),type(indHFO));
            end
            obj.plotThisTimeStamp(obj.nowCursor);
        end
        %%

        function delete(obj)
            type    = obj.newEventType;
            label   = obj.streamHandle.event.label;
            % indices = cell2mat(strfind(label,type));
            indices = strncmp(label, type, 3);
            % indices(indices==0) = [];
            latency = obj.streamHandle.event.latencyInFrame(indices);
            type = obj.streamHandle.event.label(indices);
            if isempty(latency), return;end
            
            channel = zeros(length(latency),1);
            for it=1:length(latency)
                ind = find(type{it} == '_');
                channel(it) = str2double(type{it}(ind(end)+1:end));
            end
                        
            EEG = evalin('base','EEG');

            I = ismember({EEG.event.type},obj.newEventType);
            if sum(I) == length(latency) && all(cell2mat({EEG.event(I).latency}) == latency), return;end
            EEG = pop_editeventvals(EEG,'delete',find(I));
            
            disp('Saving new events...');
            EEG = eeg_addnewevents(EEG, {latency}, {obj.newEventType},{'channel'},{channel});
            assignin('base','EEG',EEG);
            eeglab redraw
            disp('Done!!!')
        end
        %%
        
        function obj = changeSettings(obj)
            
            % store the current latency
            currentLatency = obj.nowCursor;
            
            sg = sign(double(obj.speed<1));
            sg(sg==0) = -1;
            speed1 = [num2str(-sg*obj.speed^(-sg)) 'x'];
            speed2 = [1/5 1/4 1/3 1/2 1 2 3 4 5];
            
            prefObj = [...
                PropertyGridField('gain',obj.gain,'DisplayName','Channel gain','Description','')...
                PropertyGridField('step',obj.step,'DisplayName','Page scroll by << >> in sec','Description','Page Scroll << >> in second.')...
                PropertyGridField('channels',1:obj.streamHandle.numberOfChannels,'DisplayName','Channels to plot','Description','This field accept matlab code returning a subset of channels, for instance use: ''setdiff(1:10,[3 5 7])'' to plot channels from 1 to 10 excepting 3, 5, and 7.')...
                PropertyGridField('speed',speed1,'Type',PropertyType('char','row',{'-5x','-4x','-3x','-2x','1x','2x','3x','4x','5x'}),'DisplayName','Speed','Description','Speed of the play mode.')...
                PropertyGridField('windowWidth',obj.windowWidth,'DisplayName','Window width','Description','Size of the segment to plot in seconds.')...
                PropertyGridField('normalizeFlag',obj.normalizeFlag,'DisplayName','Normalize channels','Description','Divides each channels by its standard deviation within the segment to plot.')...
                PropertyGridField('showChannelNumber',obj.showChannelNumber,'DisplayName','Show channel number or label','Description','')...
                PropertyGridField('onscreenDisplay',obj.onscreenDisplay,'DisplayName','Show events','Description','')...
                PropertyGridField('newEventType',obj.newEventType,'DisplayName','Type of new events','Description','')...
                PropertyGridField('colormap',obj.colormap,'Type',PropertyType('char','row',{'lines','eegplot'}),'DisplayName','Colormap','Description','')...
                PropertyGridField('markerType',obj.markerType,                                  'DisplayName','marker type', 'Type',PropertyType('char','row',{'+', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', 'hexagram'}),'Description','Select event marker type.')...
                PropertyGridField('markerSize',obj.markerSize,                                  'DisplayName','marker size',                  'Description','Enter the size of the event marker.')...
                PropertyGridField('markerEdgeLineWidth',obj.markerEdgeLineWidth,                'DisplayName','marker edge line width',       'Description','Enter the edge line width of the event marker.')...
                PropertyGridField('markerEdgeColor', obj.markerEdgeColor,                       'DisplayName','marker edge color',            'Description','Enter the edge color of the event marker.')...
                PropertyGridField('cursorType',obj.cursorType,                                  'DisplayName','cursor type', 'Type',PropertyType('char','row',{'+', 'o', '*', '.', 'x', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', 'hexagram'}),'Description','Select cursor type.')...
                PropertyGridField('cursorSize',obj.cursorSize,                                  'DisplayName','cursor size',                    'Description','Enter the size of the cursor.')...
                PropertyGridField('cursorLineWidth',obj.cursorLineWidth,                        'DisplayName','cursor line width',              'Description','Enter the line width of the cursor.')...
                PropertyGridField('markerSelectErrorMarginRate',obj.markerSelectErrorMarginRate,'DisplayName','error margin in marker deletion','Description','Enter the error margin [% of page width] when selecting markers for deletion.')...
                PropertyGridField('mouseWheelMoveInSec',obj.mouseWheelMoveInSec,                'DisplayName','Page scroll by mouse wheel in sec', 'Description','')...
                ];
                % PropertyGridField('labels', obj.streamHandle.event.uniqueLabel,'DisplayName','Show only a subset of events','Description','')...
                
            % create figure
            f = figure('MenuBar','none','Name','Preferences','NumberTitle', 'off','Toolbar', 'none');
            position = get(f,'position');
            set(f,'position',[position(1:2) 385 424]);
            % select the main plot
            if gca ~= obj.axesHandle
                axes(obj.axesHandle);
            end
            g = PropertyGrid(f,'Properties', prefObj,'Position', [0 0 1 1]);
            uiwait(f); % wait for figure to close
            val = g.GetPropertyValues();
            
%             obj.eventObj = event;
%             if isfield(val,'labels')
%                 if ~isempty(val.labels)
%                     for it=1:length(val.labels)
%                         latency = obj.streamHandle.event.getLatencyForEventLabel(val.labels{it});
%                         obj.eventObj = obj.eventObj.addEvent(latency,val.labels{it});
%                     end
%                 end
%             end

            obj.gain                = val.gain;
            obj.step                = val.step;
            obj.speed               = speed2(ismember({'-5x','-4x','-3x','-2x','1x','2x','3x','4x','5x'},val.speed));
            obj.windowWidth         = val.windowWidth;
            obj.normalizeFlag       = val.normalizeFlag;
            obj.showChannelNumber   = val.showChannelNumber;
            obj.channelIndex        = val.channels;
            obj.onscreenDisplay     = val.onscreenDisplay;
            obj.newEventType        = val.newEventType;
            obj.changeColormap(val.colormap);
            obj.markerType          = val.markerType;
            obj.markerSize          = val.markerSize;
            obj.markerEdgeLineWidth = val.markerEdgeLineWidth;
            obj.markerEdgeColor     = val.markerEdgeColor;
            obj.cursorType          = val.cursorType;
            obj.cursorSize          = val.cursorSize;
            obj.cursorLineWidth     = val.cursorLineWidth;
            obj.markerSelectErrorMarginRate = val.markerSelectErrorMarginRate;
            obj.mouseWheelMoveInSec = val.mouseWheelMoveInSec;
            
            figure(obj.figureHandle);
            
            % resume from the current latency
            obj.createGraphicObjects(currentLatency);
        end
        %%
        
        function plotThisTimeStamp(obj,nowCursor)
            obj.onscreenDisplay = 0;
            
            % draw line plots
            plotThisTimeStamp@eegplotNGHandle(obj,nowCursor);
            
            % update scale
            scaleHandle        = findall(obj.figureHandle, 'Tag', 'axesForScale');          % 10/7/2015 Makoto
            scaleDisplayHandle = findall(obj.figureHandle, 'Tag', 'scaleAmplitudeDisplay'); % 10/7/2015 Makoto
            set(scaleHandle,        'Units', 'normalized');
            set(scaleDisplayHandle, 'Units', 'normalized');
            axes(scaleHandle);
            cla
            mainAxesPosition = get(obj.axesHandle, 'Position');
            mainAxesOuterPosition = get(obj.axesHandle, 'OuterPosition');
%             scaleAxesPosition = get(scaleHandle, 'Position');
%             scaleAxesOuterPosition = get(scaleHandle, 'OuterPosition');
            axes(scaleHandle);
            ylim([0 1])
            numChans = get(obj.axesHandle, 'YTick');
            line([0 1], [0.0001 0.0001], 'Color', [0 0 0]);
            line([0 1], [1/(length(numChans)+1) 1/(length(numChans)+1)], 'Color', [0 0 0]);
            line([0.5 0.5], [0.001 1/(length(numChans)+1)], 'Color', [0 0 0]);
            
            set(scaleHandle, 'Position',...
                [mainAxesPosition(2)+mainAxesPosition(4)-0.015... % left end
                mainAxesPosition(2)...                            % bottom end
                0.02...                                           % width
                mainAxesPosition(4)...                            % height
                ], 'Color', [0.93 0.96 1], 'XTick', [], 'YTick', []);
            
            set(scaleDisplayHandle, 'Position',...
                [mainAxesPosition(2)+mainAxesPosition(4)-0.015+0.02... % left end
                mainAxesPosition(2)...                            % bottom end
                0.03...                                           % width
                0.015...                                           % height
                ], 'backgroundColor', [1 1 1]);
            
            if length(numChans)>1
                set(scaleDisplayHandle, 'string', sprintf('%d', round(diff(numChans([1 2])))));
            else
                set(scaleDisplayHandle, 'string', sprintf('%s', '(N.A.)'));
%                 [~,t1] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor-obj.windowWidth/2)));
%                 [~,t2] = min(abs(obj.streamHandle.timeStamp(obj.timeIndex) - (obj.nowCursor+obj.windowWidth/2)));
%                 rawdata = obj.streamHandle.mmfObj.Data.x(obj.timeIndex(t1:t2),obj.channelIndex);
%                 set(scaleDisplayHandle, 'string', sprintf('%d',round(std(rawdata)*6))); % 
            end
                
            
            
            % return if no event markers to plot
            if isempty(obj.newEventType)
                return
            end
            
            % get all event markers stored
            allEventLabels         = obj.streamHandle.event.label;
            allEventLatencyInFrame = obj.streamHandle.event.latencyInFrame;
            
            % check which ones are HFO event markers
            hfoEventIdx = strncmp(allEventLabels,obj.newEventType,3);
            
            % return if no HFO event markers
            if ~any(hfoEventIdx)
                return
            end
            
            % delete non-HFO event info
            allEventLabels(~hfoEventIdx)         = [];
            allEventLatencyInFrame(~hfoEventIdx) = [];
            
            % get channel index
            channelIdx = zeros(sum(hfoEventIdx),1);
            for n = 1:length(channelIdx)
                underscorePositions = find(allEventLabels{n}=='_');
                channelIdx(n) = str2double(allEventLabels{n}(underscorePositions(end)+1:end)); 
            end
            
            % select markers that are included in the current plotting window       
            plotWindowLatencyRange  = [obj.nowCursor-obj.windowWidth/2 obj.nowCursor+obj.windowWidth/2];
            allEventLatencyInSecond = obj.streamHandle.timeStamp(allEventLatencyInFrame);
            goodEventLatencyCriterion = find(allEventLatencyInSecond > plotWindowLatencyRange(1) & allEventLatencyInSecond < plotWindowLatencyRange(2));
            goodEventChannelCriterion = find(ismember(channelIdx, obj.channelIndex))';
            goodEventFinalIdx         = intersect(goodEventLatencyCriterion,goodEventChannelCriterion);
            
            % return if no markers left
            if ~any(goodEventFinalIdx)
                return
            end
            
            % extract values for the finally selected event markers
            goodEventLabels              = allEventLabels(goodEventFinalIdx);
            goodEventLatencyInFrameList  = allEventLatencyInFrame(goodEventFinalIdx);
            goodEventLatencyInSecondList = allEventLatencyInSecond(goodEventFinalIdx);
            goodEventChannelList         = channelIdx(goodEventFinalIdx)';
            [~,objOrder] = ismember(goodEventChannelList,obj.channelIndex);
            goodEventColorList           = obj.colorInCell(objOrder);
            tmpAmpList = get(obj.gObjHandle(objOrder),'YData');
            if iscell(tmpAmpList) % if multiple event markers
                tmpAmpList = cell2mat(tmpAmpList);
            end
            
            % obtain the current window start latency in frame
            currentWindowLatencyInSec = get(obj.gObjHandle(objOrder),'XData');
            if iscell(currentWindowLatencyInSec) % if multiple event markers
                currentWindowLatencyInSec = cell2mat(currentWindowLatencyInSec);
            end
            subtractResult = bsxfun(@minus, currentWindowLatencyInSec, goodEventLatencyInSecondList');
            [~,goodEventScreenFrame] = min(abs(subtractResult),[],2);
            goodEventChannelAmp = zeros(size(goodEventScreenFrame));
            for n = 1:length(goodEventChannelAmp)
                goodEventChannelAmp(n) = tmpAmpList(n,goodEventScreenFrame(n));
            end
            
            % delete exisitng markers to refresh
            try delete(obj.eventMarker);end %#ok
            
            % drawing markers
            hold(obj.axesHandle,'on');
            obj.eventMarker = zeros(length(goodEventFinalIdx),1);
            for n = 1:length(obj.eventMarker)
                obj.eventMarker(n) = obj.plotEventMarker(goodEventLatencyInSecondList(n),...
                    goodEventChannelAmp(n), goodEventLabels(n), goodEventColorList{n}, goodEventFinalIdx(n));
            end
            hold(obj.axesHandle,'off');
        end
        %%
        
        function markerHandle = plotEventMarker(obj,timeStamp,value,eventLabel, lineColor, latestChannel)
            markerHandle = plot(timeStamp, value, obj.markerType, 'MarkerSize', obj.markerSize, 'linewidth', obj.markerEdgeLineWidth, 'MarkerEdgeColor', obj.markerEdgeColor, 'MarkerFaceColor', lineColor, 'tag', ['diamond_ch' num2str(latestChannel)],'userdata',eventLabel,'Parent',obj.axesHandle, 'HitTest', 'off');
        end
        %%
        
        function onMouseWheelMove(obj,~,eventObj)
            % step = -10*obj.step*(eventObj.VerticalScrollCount*eventObj.VerticalScrollAmount)/obj.streamHandle.samplingRate;
            step = obj.mouseWheelMoveInSec*sign(eventObj.VerticalScrollCount*eventObj.VerticalScrollAmount);
            plotStep(obj,step);
        end
     end
end