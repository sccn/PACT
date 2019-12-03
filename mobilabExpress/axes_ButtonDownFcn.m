% axes_ButtonDownFcn()

% History
% 03/29/2014 Makoto. Debug. Now uses browserObj.channelIndex(apparentChannelOrder).
% 12/21/2012 Alejandro and Makoto. Created.

% Copyright (C) 2012 Ajelandro Ojeda, Makoto Miyakoshi, JSPS/SCCN,INC,UCSD;
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

function axes_ButtonDownFcn(hObject,~,~)
browserObj = get(get(get(hObject,'parent'),'parent'),'userData');
if isempty(browserObj)
    browserObj = get(get(hObject,'parent'),'userData');
end

cursorX = get(browserObj.cursorHandle.eventCursor, 'XData');
if cursorX < 0
    cursorX = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% determine error margin %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xMaxError = max(get(browserObj.axesHandle, 'XLim') * browserObj.markerSelectErrorMarginRate * 1/100);

switch get(browserObj.figureHandle, 'selectiontype')
%%%%%%%%%%%%%%%%%%%%%
%%% if left click %%%
%%%%%%%%%%%%%%%%%%%%%
    case 'normal'
        %- insert event
        apparentChannelOrder = get(browserObj.sliderHandle,'userdata');
        userdata = get(browserObj.axesHandle,'userdata');
        userdata = [userdata browserObj.channelIndex(apparentChannelOrder)];
        set(browserObj.axesHandle,'userdata',userdata);

        if ismac
            [~,uuid] = system('uuidgen');
        else
            uuid = java.util.UUID.randomUUID;
        end
        
        uuid = char(uuid);
        latency = cursorX;
        [~,latency] = min(abs(browserObj.streamHandle.timeStamp-latency));
        type = [browserObj.newEventType '_' uuid '_' num2str(browserObj.channelIndex(apparentChannelOrder))];

        browserObj.streamHandle.event = browserObj.streamHandle.event.addEvent(latency,type);
        %-
        voltage = get(browserObj.gObjHandle(apparentChannelOrder),'ydata');
        time = get(browserObj.gObjHandle(apparentChannelOrder),'xdata');
        [~,loc] = min(abs(time - cursorX));

        lineColor = get(browserObj.gObjHandle(apparentChannelOrder),'Color');

        hold(browserObj.axesHandle,'on');
        browserObj.eventMarker(end+1) = browserObj.plotEventMarker(time(loc),voltage(loc),type, lineColor, apparentChannelOrder);
        hold(browserObj.axesHandle,'off');
        
%%%%%%%%%%%%%%%%%%%%%%
%%% if right click %%%
%%%%%%%%%%%%%%%%%%%%%%
    case 'alt' % left mouse button click
    
        apparentChannelOrder = get(browserObj.sliderHandle,'userdata');
        dobj = findobj(browserObj.figureHandle,'tag',['diamond_ch' num2str(apparentChannelOrder)]);
        if ~isempty(dobj)
            if iscell(get(dobj,'XData'))
                xDistance = abs(cell2mat(get(dobj,'XData')) - cursorX);
            else
                xDistance = abs(get(dobj,'XData') - cursorX);
            end

            xError    = min(xDistance);
            if xError <= xMaxError
                [~,dobjIndex] = min(xDistance);
                type = get(dobj(dobjIndex),'userdata');
                delete(dobj(dobjIndex));
                browserObj.streamHandle.event = browserObj.streamHandle.event.deleteAllEventsWithThisLabel(type);   
                disp(['Event deleted: Distance to the closest marker ' num2str(xError) ' < error margin ' num2str(xMaxError)])
            else
                disp(['Event not deleted: Distance to the closest marker ' num2str(xError) ' > error margin ' num2str(xMaxError)])
            end
        else
            disp(['No event marker found in ch' num2str(browserObj.channelIndex(apparentChannelOrder))])
        end
end
end