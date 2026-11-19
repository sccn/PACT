function test_pact_circstat
%TEST_PACT_CIRCSTAT Check the bundled CircStat functions that PACT uses.
% CircStat (P. Berens, J Stat Softw 31(10), 2009) is bundled unmodified in
% external/CircStat. This test confirms that the bundled copy is the one on
% the path and that its results match closed-form definitions for the
% column-vector inputs PACT passes. No EEGLAB is required.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
bundled=fullfile(root,'external','CircStat');
for name={'circ_r','circ_mean','circ_rad2ang','circ_rtest','circ_plot','circ_std'}
    location=which(name{1});
    assert(startsWith(lower(location),lower(bundled)),...
        '%s resolves to %s, not the bundled CircStat.',name{1},location);
end

rng(17,'twister');
alpha=mod(.9+.6*randn(257,1),2*pi)-pi;   % column vector, as in PACT
n=numel(alpha);
z=mean(exp(1i*alpha));
r=abs(z);
assert(abs(circ_r(alpha)-r)<1e-12);
assert(abs(circ_mean(alpha)-angle(z))<1e-12);
assert(abs(circ_std(alpha)-sqrt(2*(1-r)))<1e-12);
[~,s0]=circ_std(alpha);
assert(abs(s0-sqrt(-2*log(r)))<1e-12);
assert(abs(circ_rad2ang(pi/3)-60)<1e-12);
[p,zStat]=circ_rtest(alpha);
assert(abs(zStat-n*r^2)<1e-9);
R=n*r;
expectedP=exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));   % Zar (1999), eq. 27.4
assert(abs(p-expectedP)<1e-12);

% Row and column inputs agree (upstream default dimension is the first
% non-singleton one; the 2012a release used dimension 1).
assert(abs(circ_r(alpha.')-r)<1e-12);
assert(abs(circ_mean(alpha.')-angle(z))<1e-12);

% circ_plot draws the polar histogram used by pac_vis_angHistPolar.
visibility=get(groot,'defaultFigureVisible');
restore=onCleanup(@()set(groot,'defaultFigureVisible',visibility));
set(groot,'defaultFigureVisible','off');
figureHandle=figure;
closeFigure=onCleanup(@()close(figureHandle));
handle=circ_plot(alpha,'hist',[],12,true,true,'linewidth',3,'color','r');
assert(~isempty(handle) && all(isgraphics(handle)));
fprintf('test_pact_circstat: PASS\n');
end
