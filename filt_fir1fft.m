% filt_fir1fft() - fast filter using fir1 and fftfilt with kaiserord() 
%                  FIR order estimator.
% Usage: [output, order, b] = filt_fir1fft(input, loLim, hiLim, samplingRate);
%        [~, order, b] = filt_fir1fft([], loLim, hiLim, samplingRate) returns the
%        design only. The output is zero-phase (delay compensated); its first and
%        last order/2 samples depend on the reflected padding (edge transients).
% Note:  Transition band-width is fixed to be 20% of the lo/liLim; pass- and
%        stopband deviation 0.01 (40 dB). Until November 2026 the deviation was
%        0.1: kaiserord then selects beta = 0 (a rectangular window) and the
%        measured passband loss reached 2.1 dB, violating the stated spec.

% Author: Makoto Miyakoshi JSPS/SCCN,INC,UCSD
% History:
% 6/20/2013 ver 1.0 by Makoto. Created.

% Copyright (C) 2013, Makoto Miyakoshi SCCN,INC,UCSD
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

function [output, n, b] = filt_fir1fft(input, loLim, hiLim, samplingRate)

if     loLim == 0 % low-pass filter
    f = [hiLim hiLim*1.2];
    a   = [1 0];
    dev = [0.01 0.01];
    freqStr = sprintf('LPF at %2.2fHz', hiLim);
elseif hiLim == 0 % high-pass filter
    f = [loLim*0.8 loLim];
    a   = [0 1];
    dev = [0.01 0.01];
    freqStr = sprintf('HPF at %2.2fHz', loLim);
else              % band-pass filter
    f = [loLim*0.8 loLim hiLim hiLim*1.2];
    a   = [0 1 0];
    dev = [0.01 0.01 0.01];
    freqStr = sprintf('BPF between %2.2f-%2.2fHz', loLim, hiLim);
end
designOnly = isempty(input) && nargout > 1;
if ~ismatrix(input) || ~isnumeric(input) || (isempty(input) && ~designOnly) || ...
        ~isscalar(samplingRate) || ~isfinite(samplingRate) || samplingRate<=0 || ...
        ~isscalar(loLim) || ~isscalar(hiLim) || loLim<0 || hiLim<0 || ...
        (loLim==0 && hiLim==0)
    error('Invalid FIR filter input or parameters.');
end
if any(f <= 0) || any(f >= samplingRate/2) || any(diff(f)<=0)
    error('FIR filter transition edges must be ordered and within Nyquist.');
end
[n,Wn,beta,ftype] = kaiserord(f,a,dev,samplingRate);
% Use the Kaiser window requested by kaiserord, not fir1's Hamming default.
% Even FIR order ensures an integral group delay.
n = 2*ceil(n/2);
% kaiserord's order is an estimate and can fall short when transition and
% pass bands are of similar width. Verify the response and lengthen the
% filter until the deviation spec holds.
[n,b] = verifiedDesign(n,Wn,beta,ftype,f,a,dev(1),samplingRate);
if designOnly
    output = [];
    return
end
nSamples = size(input,1);
if n+1 >= nSamples
    error('Filter order %d is too long for %d data samples.',n,nSamples);
end
% Boundary reflection about each endpoint, preserving its local baseline.
prefix = 2*input(1,:) - input(n+1:-1:2,:);
suffix = 2*input(end,:) - input(end-1:-1:end-n,:);
padded = [prefix; input; suffix];
filtered = fftfilt(b,padded);
first = n+n/2+1;
output = filtered(first:first+nSamples-1,:);
fprintf('filt_fir1fft: %s, order %d, TBW 20%%, PB/SB deviation 0.01.\n', freqStr, n)

function [n,b] = verifiedDesign(n,Wn,beta,ftype,f,a,dev,fs)
start = n;
while true
    b = fir1(n,Wn,ftype,kaiser(n+1,beta));
    [H,freq] = freqz(b,1,2^16,fs);
    gain = abs(H);
    if numel(f) == 2
        if a(1) == 1, pass = freq<=f(1); stop = freq>=f(2);
        else,         pass = freq>=f(2); stop = freq<=f(1); end
    else
        pass = freq>=f(2) & freq<=f(3); stop = freq<=f(1) | freq>=f(4);
    end
    if max(abs(gain(pass)-1)) <= dev && max(gain(stop)) <= dev
        return
    end
    if n > 4*start
        error('filt_fir1fft: could not meet the %.3g deviation spec (order %d).',dev,n);
    end
    n = n + 2*max(1,round(.02*n/2));
end
