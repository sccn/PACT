% onMouseMoveCallBack()

% History
% 10/08/2015 Makoto. Amplitude Scale added.
% 03/29/2014 Makoto. Debug. [~,cursor_y]=ismember(cursor_y,browserObj.channelIndex)
% 03/22/2013 Alejandro and Makoto. Created.

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

function onMouseMoveCallBack(hObject,~,~)

browserObj = get(hObject,'userData');
pos = get(browserObj.axesHandle, 'CurrentPoint');

x_lims = get(browserObj.axesHandle, 'XLim');
y_lims = get(browserObj.axesHandle, 'YLim');

x = pos(1,1);   % x position of mouse location
y = y_lims(2) - pos(1,2);    % y position of mouse location

if x >= x_lims(1) && x <= x_lims(2)  && y >= y_lims(1) && y <= y_lims(2)
    
    channels = 1:length(browserObj.channelIndex);
    nch = length(channels);
    [~,loc] = min(abs(channels -  y/(y_lims(2)/nch)));
    cursor_y = browserObj.channelIndex(channels(loc));
    
    [~,cursor_y]=ismember(cursor_y,browserObj.channelIndex); % Makoto 03/28/2014
    
    [~,t] = min(abs(browserObj.streamHandle.timeStamp - x));
    val = browserObj.streamHandle.mmfObj.Data.x(t,cursor_y);
    yData = get((browserObj.gObjHandle(cursor_y)),'ydata');
    xData = get((browserObj.gObjHandle(cursor_y)),'xdata');
    [~,loc] = min(abs(xData - x));
    pointOnChannel = yData(loc);
    
    set(browserObj.sliderHandle,'userdata',cursor_y);
    
    try delete(browserObj.cursorHandle.eventCursor);end %#ok
    hold(browserObj.axesHandle,'on');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% re-select the main plot to make the amplitude scale work %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if gca ~= browserObj.axesHandle
        axes(browserObj.axesHandle);
    end
    
    browserObj.cursorHandle.eventCursor = plot(x, pointOnChannel, browserObj.cursorType,...
        'linewidth', browserObj.cursorLineWidth, 'MarkerSize', browserObj.cursorSize,...
        'MarkerEdgeColor', get(browserObj.gObjHandle(cursor_y),'Color'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% HitTest off to make non-targeted objects unclickable %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(browserObj.cursorHandle.eventCursor,'HitTest','off');
    set(browserObj.gObjHandle,{'LineWidth'},{0.5}, {'HitTest'}, {'off'});
    set(browserObj.gObjHandle(cursor_y),'LineWidth',1.5,'HitTest','off'); % This changes the line width of the currently selected channel
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Make the entire axis clickable %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(browserObj.axesHandle,'ButtonDownFcn',@axes_ButtonDownFcn);
    hold(browserObj.axesHandle,'off');
end