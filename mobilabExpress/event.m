% Class that defines the behavior and functionality of MoBI events
%
% Author: Nima Bigdely Shamlo, Alejandro Ojeda, SCCN, INC, UCSD
classdef event % < hedManager
    properties
        latencyInFrame = []; %  1 x N vector containing latencies for N events
        label = {};  %  1 x N vector containing labels (strings describing event types) for N events.
        hedTag;
    end
    properties(Dependent)
        uniqueLabel;
    end
    properties (Hidden = true)
        value = [];
    end
    methods
        %%
        function obj = event(eventChannel,label)
            if nargin < 1, eventChannel = 0;end
            if isstruct(eventChannel) || isa(eventChannel,'event')
                obj.latencyInFrame = eventChannel.latencyInFrame;
                if isfield(eventChannel,'hedTag')
                     obj.hedTag = eventChannel.hedTag;
                else obj.hedTag = eventChannel.label;
                end
                obj.label = eventChannel.label;
                if isfield(eventChannel,'value'), obj.value = eventChannel.value;end
            else 
                obj.latencyInFrame = find(eventChannel);
                N = length(obj.latencyInFrame);
                if nargin < 2
                     obj.hedTag = mat2cell(num2str(eventChannel(obj.latencyInFrame)),ones(N,1));
                else obj.hedTag = label;
                end
                indBoundary = isnan(eventChannel);
                [~,ind] = ismember(obj.latencyInFrame,find(indBoundary));
                obj.hedTag(ind) = {'boundary'};
                obj.label = obj.hedTag; 
            end
        end
        function uniqueLabel = get.uniqueLabel(obj), uniqueLabel = unique(obj.label);end
        function obj = addEventFromChannel(obj,eventChannel,offset,label)
            
            if nargin < 3, offset = 0;end
            eventChannel = eventChannel(:);
            tmp = find(eventChannel);
            N = length(tmp);
            obj.latencyInFrame(end+1:end+N) = offset+tmp;
            if nargin < 4, 
                if isempty(obj.hedTag)
                    obj.hedTag = cell(N,1);
                    for it=1:N
                        obj.hedTag{it} = strtrim(num2str(eventChannel(tmp(it))));
                        obj.label{it} = obj.hedTag{it};
                    end
                else
                    for it=1:N
                        obj.hedTag{end+1} = strtrim(num2str(eventChannel(tmp(it))));
                        obj.label{end+1} = strtrim(num2str(eventChannel(tmp(it))));
                    end
                end
            else
                obj.hedTag(end+1:end+N) = label;
                obj.label(end+1:end+N) = label;
            end
            indBoundary = isnan(eventChannel);
            [~,ind] = ismember(obj.latencyInFrame,offset+find(indBoundary));
            obj.hedTag(logical(ind)) = {'boundary'}; 
            obj.label(logical(ind)) = {'boundary'}; 
        end
        %%
        function obj = addEvent(obj,latencyInFrame,label,varargin)
            if nargin < 3, error('prog:input','Not enough input arguments.');end
            N = length(latencyInFrame);
            obj.latencyInFrame(end+1:end+N) = latencyInFrame;
            shortLabel = label;

            if iscell(label)
                obj.hedTag(end+1:end+N) = label;
                if event.isHedTag(label{1});
                    for it=1:length(label)
                        loc = find(label{it} == '/');
                        if length(loc) > 2, loc = loc(2);
                        else loc = length(label{it})+1;
                        end
                        shortLabel{it} = label{it}(1:loc-1);
                    end
                end
                obj.label(end+1:end+N) = shortLabel;
            else
                obj.hedTag(end+1:end+N) = {label};
                if event.isHedTag(label);
                    loc = find(label == '/');
                    if length(loc) > 2, loc = loc(2);
                    else loc = length(label)+1;
                    end
                    obj.label(end+1:end+N) = {label(1:loc(end)-1)};
                else obj.label(end+1:end+N) = {label};
                end 
            end
            
            if nargin > 3
                name = varargin(1:2:end); 
                val = varargin(2:2:end);%#ok
                if isempty(obj.value)
                    for it=1:length(name), eval(['obj.value(1).' name{it} '= val{it}(:);']);end
                    obj.value(1).label = shortLabel; 
                else
                    ind = ismember({obj.value.label},shortLabel);
                    if ~any(ind)
                        n = length(obj.value);
                        for it=1:length(name), eval(['obj.value(n+1).' name{it} '= val{it}(:);']);end
                        obj.value(n+1).label = shortLabel;
                    else
                        ind = find(ind);
                        ind = ind(1);%#ok
                        for it=1:length(name), eval(['obj.value(ind).' name{it} '(end+1:end+N)= val{it}(:);']);end
                    end
                end
            end
        end
        %%
        function obj = eeglab2event(obj,EEG)
            obj.latencyInFrame = {EEG.event.latency};
            obj.hedTag = {EEG.event.type};
            % eventType must contain strings
            for it = 1:length(obj.latencyInFrame)
                if isnumeric(obj.hedTag{it}), obj.hedTag{it} = num2str(obj.hedTag{it});end                
            end
        end
        %%
        function EEG = event2eeglab(obj,EEG) 
            N = length(obj.label);
            if N
                if N == 1
                    latencies = {obj.latencyInFrame};
                elseif N == 2
                    latencies = {obj.latencyInFrame(1),obj.latencyInFrame(2)};
                else
                    latencies = num2cell(obj.latencyInFrame,[N 1]);
                end
                EEG = eeg_addnewevents(EEG, latencies, obj.label,{'hedTag'},obj.hedTag);
            end
        end
        %%
        function metadata = saveobj(obj)
            metadata.latencyInFrame = obj.latencyInFrame;
            metadata.label = obj.label;
            metadata.value = obj.value;
            metadata.hedTag = obj.hedTag;
        end
        %%
        function numberOfOccurancesForEachEvent = getNumberOfOccurancesForEachEvent(obj)            
            N = length(obj.uniqueLabel);
            numberOfOccurancesForEachEvent = zeros(1,N);
            for it=1:N, numberOfOccurancesForEachEvent(it) = sum(ismember(obj.label,obj.uniqueLabel{it}));end            
        end
        %%
        function plotNumberOfOccurancesForEachEvent(obj)
            numberOfOccurancesForEachEvent = getNumberOfOccurancesForEachEvent(obj);
            figure;
            barh(numberOfOccurancesForEachEvent);
            set(gca,'ytick', 2:length(obj.uniqueLabel),  'ytickLabel', obj.uniqueLabel(2:end))
            ylabel('Events');
            xlabel('Number of occurences');
        end
        %%
        function [eventLatencyInFrame, eventLabel] = getLatencyForEventLabel(obj, eventLabel)
            if isnumeric(eventLabel), eventLabel = num2str(eventLabel);end
            [~,loc] = ismember(obj.label,eventLabel); 
            eventLatencyInFrame = obj.latencyInFrame(logical(loc));
        end
        %%
        function obj = renameLabels(obj,label,newLabel)
           loc = find(ismember(obj.label,label));
           if isempty(loc), return;end
           for it=1:length(loc), obj.label{loc(it)} = newLabel;end
        end
        %%
        function obj = deleteEvent(obj,index)
            if length(nonzeros(index)) <= length(obj.latencyInFrame) && ~isempty(nonzeros(index))
                obj.latencyInFrame(index) = [];
                obj.label(index) = [];
                obj.hedTag(index) = [];
            end
        end
        %%
        function obj = deleteAllEventsWithThisLabel(obj,label)
            if iscell(label)
                for it=1:length(label)
                    index = ismember(obj.label,label{it});
                    obj = obj.deleteEvent(index);
                end
            else
                index = ismember(obj.label,label);
                obj = obj.deleteEvent(index);
            end
            if ~isempty(obj.value)
                ind = ismember({obj.value.label},label);
                obj.value(ind) = [];
            end
        end
    end
    methods(Static=true)
        %%
        function obj = loadobj(a)
            obj = event;
            obj.latencyInFrame = a.latencyInFrame;
            obj.label = a.label;
            obj.value = a.value;
            obj.hedTag = a.hedTag;
        end
        %%
        function val = isHedTag(strLabel)
            val = false;
            if isempty(strLabel), return;end
            %manager = hedManager;
            %val = manager.isValidHedTag(strLabel);
            if any(strLabel == '/'), val = true;end
        end
    end
    methods(Hidden=true)
        function eventId = getIdForEventLabel(obj, eventLabel), eventId = find(strcmp(obj.uniqueLabel, eventLabel));end
        %%
        function obj = interpEvent(obj,x,xi)
            yi = zeros(size(obj.latencyInFrame));
            if ~isempty(obj.latencyInFrame)
                for it=1:length(obj.latencyInFrame)
                    [~,loc] = min(sqrt((x(obj.latencyInFrame(it))-xi).^2));
                    yi(it) = loc;
                end
                obj.latencyInFrame = yi;
            end
        end
    end
end