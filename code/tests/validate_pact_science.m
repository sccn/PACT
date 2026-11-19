function validate_pact_science
% Each channel is an independent realization. Reuse noise across strengths.
rng(521,'twister');
fs=100;
n=9000;
batchSize=160;
nBatches=5;
rates=[1 2 5 20 100];
strengths=[0 .12 .5 0];
rejections=zeros(numel(strengths),numel(rates));
missing=zeros(size(rejections));
t=(0:n-1)/fs;
e=eeg_emptyset;
e.srate=fs; e.pnts=n; e.nbchan=batchSize; e.trials=1; e.xmax=t(end);
for batch=1:nBatches
    lfphase=2*pi*3*t+2*pi*rand(batchSize,1)+cumsum(.07*randn(batchSize,n),2);
    lf=cos(lfphase);
    amp=exp(.45*filter(sqrt(1-.96^2),[1 -.96],randn(batchSize,n),[],2));
    carrier=cos(2*pi*28*t+2*pi*rand(batchSize,1));
    noise=.15*randn(batchSize,n);
    for scenario=1:numel(strengths)
        e.data=lf+amp.*(1+strengths(scenario)*lf).*carrier+noise;
        if scenario==4
            e.data(:,n/2:end)=1.7*e.data(:,n/2:end);
            for ch=1:batchSize
                ix=randi([500 n-500],1,8);
                e.data(ch,ix)=e.data(ch,ix)+8*randn(1,8);
            end
        end
        for k=1:numel(rates)
            rate=rates(k);
            evalc('o=pac_man(e,''lfoPhase'',[2 4],''hfoAmp'',[20 36],''hfoTopRatio'',rate,''numSurro'',199)');
            pv=o.pac.uncorrected.MIpval;
            rejections(scenario,k)=rejections(scenario,k)+sum(pv<=.05);
            missing(scenario,k)=missing(scenario,k)+sum(~isfinite(pv));
        end
    end
end
for scenario=1:numel(strengths)
    for k=1:numel(rates)
        fprintf('RESULT %d %g %d/%d missing %d\n',scenario,rates(k),...
            rejections(scenario,k),batchSize*nBatches,missing(scenario,k));
    end
end

% Conditional rotation null: block weights are fixed and observed block
% orientations are iid uniform. These are the randomization assumptions.
rng(290,'twister');
blockCount=20;
blockLength=100;
indices=1:blockLength:blockCount*blockLength;
weights=exp(linspace(-2,2,blockCount));
valid=true(1,blockCount*blockLength);
nullRejections=0;
for replicate=1:1000
    phase=zeros(size(valid));
    phase(indices)=2*pi*rand(1,blockCount);
    mi=abs(mean(weights.*exp(1i*phase(indices))));
    null=pac_shiftSurrogate(phase,weights,indices,valid,199,100,5);
    nullRejections=nullRejections+((1+sum(null>=mi))/200<=.05);
end
fprintf('ROTATION NULL %d/1000\n',nullRejections);

% Phase-locked windows can concentrate phases with an unmodulated carrier.
e.nbchan=1;
e.data=cos(2*pi*4*t)+cos(2*pi*28*t);
centers=501:500:n-500;
e.event=struct('type','stim','latency',num2cell(centers));
opts={'lfoPhase',[3 5],'hfoAmp',[20 36],'hfoTopRatio',100,...
    'hfoPool',4,'windowLength',40,'numSurro',199};
evalc('o=pac_man(e,opts{:});');
fprintf('PHASE-LOCKED MARKER, UNMODULATED CARRIER MI %.4g p %.4g\n',...
    o.pac.mi,o.pac.uncorrected.MIpval);
end
