function dismisser = pact_dialogDismisser(deadlineSeconds)
%PACT_DIALOGDISMISSER Cancel real modal dialogs in GUI smoke tests, bounded.
%   d = pact_dialogDismisser(deadlineSeconds);   % arm before opening dialog
%   ... code that opens a real dialog and waits (inputgui, inputdlg, ...)
%   log = pact_dialogDismisserStop(d);           % disarm; returns the log
%
% A timer polls every 0.2 s for new figures that contain uicontrols. Once a
% dialog waits for the user (WaitStatus 'waiting', or visible for 1 s for
% waitfor-based dialogs such as inputgui), it presses its Cancel button
% through the button's own callback,
% matching Tag or String 'cancel' case-insensitively (MATLAB's inputdlg uses
% 'Cancel'; EEGLAB's inputgui uses 'cancel'). If a dialog has no Cancel button,
% or is still open deadlineSeconds after it was first seen, the figure is
% deleted and the event is logged as 'timeout', which unblocks uiwait and
% waitfor. Callers must treat any 'timeout' as a test failure.
%
% The timer body never throws. An exception in a TimerFcn would stop the
% timer, which is how the original harness could leave a dialog waiting forever.
%
% Log entries: struct with fields name, action ('canceled'|'timeout'),
% seconds (time from first sighting), detail.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

validateattributes(deadlineSeconds,{'numeric'},{'scalar','positive','finite'});
state.existing=findall(groot,'Type','figure');
state.deadline=deadlineSeconds;
state.clock=tic;
state.seen=struct('figure',{},'time',{},'pressed',{},'visibleSince',{});
state.log=struct('name',{},'action',{},'seconds',{},'detail',{});
dismisser=timer('ExecutionMode','fixedSpacing','Period',.2,'BusyMode','drop',...
    'Name','PACT dialog dismisser','UserData',state,'TimerFcn',@tick);
start(dismisser);
end

function tick(timerObject,~)
try
    state=timerObject.UserData;
    figures=findall(groot,'Type','figure');
    for f=reshape(figures,1,[])
        if ~isvalid(f) || any(state.existing==f), continue, end
        controls=findall(f,'Type','uicontrol');
        if isempty(controls), continue, end          % a plot, not a dialog
        k=find(arrayfun(@(s)isequal(s.figure,f),state.seen),1);
        if isempty(k)
            state.seen(end+1)=struct('figure',f,'time',toc(state.clock),'pressed',false,'visibleSince',NaN);
            k=numel(state.seen);
        end
        elapsed=toc(state.clock)-state.seen(k).time;
        if isnan(state.seen(k).visibleSince) && strcmp(get(f,'Visible'),'on')
            state.seen(k).visibleSince=toc(state.clock);
        end
        name=char(get(f,'Name'));
        button=cancelButton(controls);
        if state.seen(k).pressed && elapsed<state.deadline
            continue                                   % closing; give it time
        elseif ~isempty(button) && elapsed<state.deadline && ready(f,state.seen(k),state.clock)
            % Press only once the dialog is waiting for the user, as a user would.
            state.seen(k).pressed=true;
            try
                detail=press(button,f);
                state.log(end+1)=entry(name,'canceled',elapsed,detail);
            catch ME
                delete(f);
                state.log(end+1)=entry(name,'timeout',elapsed,['Cancel callback failed: ' ME.message]);
            end
        elseif elapsed>=state.deadline
            delete(f);
            detail='no Cancel button found';
            if ~isempty(button), detail='dialog stayed open past the deadline'; end
            state.log(end+1)=entry(name,'timeout',elapsed,detail);
        end
    end
    timerObject.UserData=state;
catch ME
    % Never let the timer die; record the failure instead.
    try
        state=timerObject.UserData;
        state.log(end+1)=entry('(timer)','timeout',NaN,ME.message);
        timerObject.UserData=state;
    catch
    end
end
end

function yes=ready(f,seen,clock)
% uiwait (inputdlg, questdlg) sets WaitStatus='waiting'. EEGLAB's inputgui
% blocks in waitfor, which sets no flag, so require 1 s of visibility.
yes=strcmp(get(f,'WaitStatus'),'waiting') || ...
    (~isnan(seen.visibleSince) && toc(clock)-seen.visibleSince>=1);
end

function button=cancelButton(controls)
button=[];
for c=reshape(controls,1,[])
    if ~strcmpi(get(c,'Style'),'pushbutton'), continue, end
    if strcmpi(strtrim(char(get(c,'Tag'))),'cancel') || ...
            strcmpi(strtrim(char(get(c,'String'))),'cancel')
        button=c;
        return
    end
end
end

function detail=press(button,figureHandle)
% Run the button's own callback as MATLAB would on a click.
detail='Cancel button callback';
callback=get(button,'Callback');
if isempty(callback)
    delete(figureHandle);
elseif ischar(callback) || isstring(callback)
    % String callbacks such as EEGLAB's 'close(gcbf)' refer to gcbf/gcbo.
    code=regexprep(char(callback),'\<gcbf\>','figureHandle');
    code=regexprep(code,'\<gcbo\>','button');
    eval(code);
elseif isa(callback,'function_handle')
    callback(button,[]);
elseif iscell(callback)
    feval(callback{1},button,[],callback{2:end});
end
% Callbacks that act on gcbf (e.g. inputdlg's delete(gcbf)) do nothing
% outside a genuine click, because gcbf cannot be set programmatically. If
% the dialog is still blocked in uiwait, close it the way a user closing the
% window would (its CloseRequestFcn); dialogs document this as Cancel.
if isvalid(figureHandle) && strcmp(get(figureHandle,'WaitStatus'),'waiting')
    close(figureHandle);
    detail='Cancel callback needs gcbf; closed via CloseRequestFcn';
end
end

function e=entry(name,action,seconds,detail)
e=struct('name',name,'action',action,'seconds',seconds,'detail',detail);
end
