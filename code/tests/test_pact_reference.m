function test_pact_reference
% Independent threshold-based multiplicity checks, including tied ranks.
rng(38); p=[rand(1,33).^3 .02 .02 0 1];
sorted=sort(p); m=numel(p);
for k=1:m
    ranks=find(sorted<=p(k)); holm(k)=min(1,max(sorted(ranks).*(m+1-ranks)));
    ranks=find(sorted>=p(k)); bh(k)=min(1,min(sorted(ranks).*m./ranks));
end
assert(max(abs(holm-pac_adjustPValues(p,'holm')))<1e-12);
assert(max(abs(bh-pac_adjustPValues(p,'fdr')))<1e-12);
n=8192; fs=128; center=n/2; x=zeros(n,1); x(center)=1;
for filter={@filt_fir1fft,@filt_fir1fft_slower}
    for limits={[0 20],[10 0],[10 20]}
        lim=limits{1}; y=filter{1}(x,lim(1),lim(2),fs);
        [~,peak]=max(abs(y)); assert(peak==center);
        assert(max(abs(y(1:end-1)-flip(y(1:end-1))))<1e-10);
        gain=abs(fft(y)); probes=gain(1+[4 15 35]*n/fs);
        expected=[lim(1)==0 true lim(2)==0];
        assert(all(abs(probes(:).'-expected)<.12));
    end
end
fprintf('test_pact_reference: PASS\n');
end
