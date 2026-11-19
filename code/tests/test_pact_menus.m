function test_pact_menus
% Exercise registered callback strings and the real plotting backends.
oldVisibility=get(groot,'defaultFigureVisible');
restoreVisibility=onCleanup(@()set(groot,'defaultFigureVisible',oldVisibility));
set(groot,'defaultFigureVisible','off');
existingFigures=findall(groot,'Type','figure');
f=figure; uimenu(f,'tag','tools'); eegplugin_pac(f,struct,struct);
menu=findall(f,'type','uimenu','label','PACT');
items=flipud(get(menu,'Children'));
assert(numel(items)==11);
rng(79); EEG=eeg_emptyset; EEG.srate=100; EEG.nbchan=2;
EEG.pnts=8000; EEG.trials=1; EEG.xmax=79.99; EEG.data=randn(2,8000);
EEG.chanlocs=struct('labels',{'A','B'});
EEG.event=struct('type','stim','latency',num2cell(500:500:7500),'chan','A');
EEG.urevent=EEG.event;
for k=1:numel(EEG.event), EEG.event(k).urevent=k; end
original=EEG;
CURRENTSET=1; ALLEEG=[EEG EEG];
% Supply deterministic answers only at the modal input boundary.
folder=tempname; mkdir(folder);
cleanup=onCleanup(@()cleanupFixtures(folder,existingFigures));
writelines(["function result=inputgui(varargin)"; ...
    "queue=getappdata(0,'PACTanswers'); result=queue{1};"; ...
    "setappdata(0,'PACTanswers',queue(2:end)); end"],fullfile(folder,'inputgui.m'));
addpath(folder,'-begin');
answers={{'[2 4]','[20 35]','20',4,1,'500','.05','99','12'},...
    {2},{'2'},{'.05',4,2},{},{},{},{},{'1'},...
    {'2 5','4','20 35','1 100','4',1,1},{'2 5','20 35','100','4',1,1}};
for k=1:numel(items)
    cb=get(items(k),'callback'); label=get(items(k),'label');
    if ~isempty(answers{k})
        before=EEG; beforeAll=ALLEEG;
        setappdata(0,'PACTanswers',{{}});
        eval(cb);
        assert(isequaln(before,EEG) && isequaln(beforeAll,ALLEEG));
        fprintf('CANCEL PASS %s\n',label);
    end
    queue={answers{k}};
    if k==1, queue{2}={}; end
    setappdata(0,'PACTanswers',queue);
    eval(cb); drawnow;
    if k==3
        assert(isequaln(ALLEEG(2).event,EEG.event));
        assert(isequaln(ALLEEG(2).urevent,EEG.urevent));
    end
    assert(isequaln(EEG.data,original.data) && isequaln(EEG.event,original.event));
    assert(isequaln(EEG.urevent,original.urevent));
    fprintf('MENU PASS %s\n',label);
end
end

function cleanupFixtures(folder,existingFigures)
rmpath(folder); clear inputgui
delete(fullfile(folder,'inputgui.m')); rmdir(folder); rmappdata(0,'PACTanswers');
delete(setdiff(findall(groot,'Type','figure'),existingFigures));
end
