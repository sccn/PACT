function validate_pact_compatibility(legacy)
% Optional legacy is a handle to the base revision's pac_man under a
% different function name. See VALIDATION.md for the extraction command.
rng(464,'twister');
e=eeg_emptyset;
e.srate=100; e.nbchan=8; e.pnts=360000; e.trials=1;
e.xmax=(e.pnts-1)/e.srate; e.data=randn(e.nbchan,e.pnts,'single');
opts={'lfoPhase',[2 4],'hfoAmp',[20 35],'hfoTopRatio',5,'hfoPool',1,...
    'whichMarker',1,'windowLength',[],'alpha',.05,'numSurro',199,'numPhaseBin',12};

start=tic;
evalc('result=pac_man(e,opts{:});');
fprintf('BENCHMARK 8 x 360000, 199 draws, %.3f seconds\n',toc(start));
assert(isequal(result.data,e.data));
assert(all(isfinite(result.pac.uncorrected.MIpval)));

if nargin==0
    return
end
e.data=double(e.data(1:2,1:12000)); e.pnts=12000; e.xmax=119.99; e.nbchan=2;
e.event=struct('type','stim','latency',500,'duration',0);
e.chanlocs=struct('labels',{'A','B'});
root=groot;
visibility=get(root,'defaultFigureVisible');
restore=onCleanup(@()set(root,'defaultFigureVisible',visibility));
set(root,'defaultFigureVisible','off');
evalc('result=pac_man(e,opts{:}); prior=legacy(e,opts{:});');
for name={'mi','vectLength','angleMean','angleStd'}
    key=name{1};
    delta=max(abs(result.pac.(key)-prior.pac.(key)));
    fprintf('LEGACY %s max absolute difference %.12g\n',key,delta);
    assert(delta<1e-12);
end
end
