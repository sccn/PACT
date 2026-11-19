function results = validate_pact_science(varargin)
%VALIDATE_PACT_SCIENCE Monte Carlo calibration and power of the PACT MI test.
%   results = validate_pact_science            % 800 channels per cell
%   results = validate_pact_science('quick')   % 160 channels per cell
%
% Requires EEGLAB (pop_eegfiltnew) and PACT on the MATLAB path. Synthetic
% continuous recordings: 90 s at 100 Hz, 160 independent channels per
% batch, fixed seeds. Each channel is one realization. The LFO is a 3 Hz
% oscillation with a random-walk phase; the "HFO" is a 28 Hz carrier with
% an autocorrelated log-normal envelope. Coupling multiplies the envelope by
% (1 + k*cos(LFO phase)). Reported rates use alpha = 0.05 on uncorrected
% p-values with 199 sign-flip draws; 95% Clopper-Pearson intervals.
%
% Scenario codes
%   N1  independent; AR(1) log-envelope, time constant 0.25 s
%   N2  independent; long bursts, time constant 1 s, heavier tails
%   N3  independent; non-sinusoidal LFO (harmonics), LFO band 2-10 Hz
%   N4  pink (1/f) Gaussian noise only, no oscillations
%   N5  N1 + broadband EMG-like bursts, slow gain drift, 1.7x gain step,
%       and 4 EEGLAB boundaries between independently generated segments
%   N6  N1, event-locked mode (hfoPool=4), events unrelated to the LFO
%   C1  N1 + 8 broadband impulses per channel (artifact confound: impulses
%       genuinely couple filtered LFO phase and HFO envelope)
%   P   N1 with coupling strength k (power)
%   PE  event-locked mode, coupling only inside +/-0.5 s event windows
%
% LEGACY columns re-test N1-N5 with the historical PACT null (random
% permutation of single phase samples), using the same HAS samples and MI.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

quick = nargin>0 && strcmpi(varargin{1},'quick');
nBatch = 5; if quick, nBatch = 1; end
fs=100; n=9000; nCh=160; nSurro=199; alpha=.05;
rates=[1 2 5 20 100];
familySize=8;
strengths=[.05 .1 .2 .3 .5];
results=struct('scenario',{},'rate',{},'rejected',{},'tested',{},...
    'missing',{},'familyRejected',{},'families',{},'legacyRejected',{},...
    'medianBlocks',{});

nulls={'N1','N2','N3','N4','N5','N6','C1'};
for s=1:numel(nulls)
    results=runScenario(results,nulls{s},0,s);
end
for k=1:numel(strengths)
    results=runScenario(results,sprintf('P k=%.2f',strengths(k)),strengths(k),100+k);
end
results=runScenario(results,'PE k=0.30',.3,200);
results=runScenario(results,'PE k=0.50',.5,201);

fprintf('\n%-11s %5s %16s %18s %8s %16s %7s\n','scenario','HAS%',...
    'rejected/tested','rate [95% CI]','missing','Holm FWER (8 ch)','blocks');
for r=results
    [rate,ci]=binofit(r.rejected,max(1,r.tested));
    line=sprintf('%-11s %5g %7d/%-8d %5.1f%% [%4.1f,%4.1f] %8d %6d/%-4d %4.1f%% %7.0f',...
        r.scenario,r.rate,r.rejected,r.tested,100*rate,100*ci(1),100*ci(2),...
        r.missing,r.familyRejected,r.families,100*r.familyRejected/max(1,r.families),...
        r.medianBlocks);
    if ~isnan(r.legacyRejected)
        [lrate,lci]=binofit(r.legacyRejected,max(1,r.tested+r.missing));
        line=[line sprintf('   LEGACY %5.1f%% [%4.1f,%4.1f]',100*lrate,100*lci(1),100*lci(2))]; %#ok<AGROW>
    end
    fprintf('RESULT %s\n',line);
end

    function results=runScenario(results,name,strength,seed)
        code=strtok(name);
        rejected=zeros(size(rates)); tested=rejected; missing=rejected;
        famRej=rejected; fams=rejected; legacy=nan(size(rates)); blocks=cell(size(rates));
        if any(strcmp(code,{'N1','N2','N3','N4','N5'})), legacy(:)=0; end
        for batch=1:nBatch
            rng(1000*seed+batch,'twister');
            [e,opts]=simulate(code,strength);
            for r=1:numel(rates)
                o=quietPac(e,[opts {'hfoTopRatio',rates(r),'numSurro',nSurro,'alpha',alpha}]);
                p=o.pac.uncorrected.MIpval;
                ok=isfinite(p);
                rejected(r)=rejected(r)+sum(p(ok)<=alpha);
                tested(r)=tested(r)+sum(ok);
                missing(r)=missing(r)+sum(~ok);
                blocks{r}=[blocks{r}; o.pac.numBlocks(ok)];
                for f=1:familySize:nCh
                    adjusted=pac_adjustPValues(p(f:f+familySize-1),'holm');
                    if any(isfinite(adjusted))
                        fams(r)=fams(r)+1;
                        famRej(r)=famRej(r)+any(adjusted<=alpha);
                    end
                end
                if ~isnan(legacy(r))
                    legacy(r)=legacy(r)+legacyRejections(o,nSurro,alpha);
                end
            end
        end
        for r=1:numel(rates)
            results(end+1)=struct('scenario',name,'rate',rates(r),...
                'rejected',rejected(r),'tested',tested(r),'missing',missing(r),...
                'familyRejected',famRej(r),'families',fams(r),...
                'legacyRejected',legacy(r),'medianBlocks',median(blocks{r})); %#ok<AGROW>
        end
        fprintf('done %s\n',name);
    end

    function [e,opts]=simulate(code,k)
        t=(0:n-1)/fs;
        e=eeg_emptyset;
        e.srate=fs; e.pnts=n; e.nbchan=nCh; e.trials=1; e.xmin=0; e.xmax=t(end);
        e.chanlocs=struct('labels',arrayfun(@(c)sprintf('S%d',c),1:nCh,'UniformOutput',false));
        opts={'lfoPhase',[2 4],'hfoAmp',[20 36],'hfoPool',1,'numPhaseBin',36};
        switch code
            case 'N2'
                e.data=oscillators(t,0,.99,.8,false);
            case 'N3'
                e.data=oscillators(t,0,.96,.45,true);
                opts{2}=[2 10];
            case 'N4'
                spectrum=fft(randn(nCh,n),[],2);
                f=[0:n/2 -(n/2-1):-1]*fs/n;
                gain=1./sqrt(max(abs(f),.5));
                x=real(ifft(bsxfun(@times,spectrum,gain),[],2));
                e.data=bsxfun(@rdivide,x,std(x,0,2));
            case 'N5'
                segments=5; len=n/segments; e.data=zeros(nCh,n);
                for g=1:segments
                    span=(g-1)*len+(1:len);
                    e.data(:,span)=oscillators(t(1:len),0,.96,.45,false);
                end
                e.event=struct('type','boundary','latency',num2cell((1:segments-1)*len+.5),'duration',0);
                drift=1+.4*sin(bsxfun(@plus,2*pi*t/37,2*pi*rand(nCh,1)));
                e.data=e.data.*drift;
                e.data(:,n/2+1:end)=1.7*e.data(:,n/2+1:end);
                for c=1:nCh
                    for b=1:6
                        len=randi([20 100]); first=randi([1 n-len]);
                        e.data(c,first:first+len-1)=e.data(c,first:first+len-1)+2*randn(1,len);
                    end
                end
            case 'C1'
                e.data=oscillators(t,0,.96,.45,false);
                for c=1:nCh
                    ix=randi([300 n-300],1,8);
                    e.data(c,ix)=e.data(c,ix)+8*randn(1,8);
                end
            case {'N6','PE'}
                centers=150:300:n-150;
                centers=centers+randi([-50 50],size(centers));
                gate=false(1,n);
                for c=centers, gate(c-50:c+50)=true; end
                e.data=oscillators(t,k*gate,.96,.45,false);
                e.event=struct('type','stim','latency',num2cell(centers),'duration',0);
                opts(5:6)={'hfoPool',4};
                opts=[opts {'whichMarker',1,'windowLength',500}];
                opts{8}=12;   % fewer phase bins: event windows hold fewer samples
            otherwise   % N1 and P
                e.data=oscillators(t,k,.96,.45,false);
        end
        e.urevent=e.event;
    end

    function x=oscillators(t,k,phi,sd,nonsinusoidal)
        m=numel(t);
        theta=bsxfun(@plus,2*pi*3*t,2*pi*rand(nCh,1))+cumsum(.05*randn(nCh,m),2);
        lf=cos(theta);
        if nonsinusoidal
            lf=lf+.6*cos(2*theta+pi/3)+.3*cos(3*theta+1);
        end
        ar=filter(sqrt(1-phi^2),[1 -phi],randn(nCh,m+2000),[],2);
        envelope=exp(sd*ar(:,2001:end));
        carrier=cos(bsxfun(@plus,2*pi*28*t,2*pi*rand(nCh,1))+cumsum(.2*randn(nCh,m),2));
        x=lf+envelope.*(1+bsxfun(@times,k,cos(theta))).*carrier+.15*randn(nCh,m);
    end
end

function o=quietPac(e,options)
% Run pac_man without filter progress output (evalc needs a non-nested scope).
evalc('o=pac_man(e,options{:});');
end

function count=legacyRejections(o,nSurro,alpha)
% Historical PACT null: phases drawn at random (without replacement) from
% all eligible samples, paired with the fixed HAS amplitudes.
eligible=find(o.pac.selectedSampleMask);
count=0;
for c=1:numel(o.pac.hfoIndex)
    index=o.pac.hfoIndex{c};
    if numel(index)<4, continue, end
    z=o.pac.analyticEEG(c,:);
    amplitude=abs(z(index));
    observed=abs(mean(z(index)));
    unitPhasor=exp(1i*angle(z(eligible)));
    null=zeros(1,nSurro);
    for s=1:nSurro
        null(s)=abs(mean(amplitude.*unitPhasor(randperm(numel(eligible),numel(index)))));
    end
    count=count+((1+sum(null>=observed))/(nSurro+1)<=alpha);
end
end
