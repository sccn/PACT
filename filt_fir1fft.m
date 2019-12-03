% filt_fir1fft() - fast filter using fir1 and fftfilt with kaiserord() 
%                  FIR order estimator.
% Usage: output = filt_fir1fft(input, loLim, hiLim, samplingRate);
% Note:  Transition band-width is fixed to be 20% of the lo/liLim, and
%        pass/stopband ripple is fixed to be 0.1.

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

function output = filt_fir1fft(input, loLim, hiLim, samplingRate)

if     loLim == 0 % low-pass filter
    f = [hiLim hiLim*1.2];
    a   = [1 0];
    dev = [0.1 0.1];
    freqStr = sprintf('LPF at %2.2fHz', hiLim);
elseif hiLim == 0 % high-pass filter
    f = [loLim*0.8 loLim];
    a   = [0 1];
    dev = [0.1 0.1];
    freqStr = sprintf('HPF at %2.2fHz', loLim);
else              % band-pass filter
    f = [loLim*0.8 loLim hiLim hiLim*1.2];
    a   = [0 1 0];
    dev = [0.1 0.1 0.1];
    freqStr = sprintf('BPF between %2.2f-%2.2fHz', loLim, hiLim);
end
[n,Wn,beta,ftype] = kaiserord(f,a,dev,samplingRate);
b = fir1(n,Wn,ftype);
% fvtool(b)
prefix  = input(n:-1:1,:)*-1;
suffix  = input(end:-1:end-n+1,:)*-1;
tmpSgnl = [prefix; input; suffix];
tmpOut1 = fftfilt(b,tmpSgnl);
tmpOut2 = tmpOut1(floor(n*1.5)+1:end,:); % align the beginning; +1 as a result of comparison with eegfilt
output  = tmpOut2(1:size(input,1),:);
fprintf('filt_fir1fft: %s, order %d, TBW 20%%, PB/SB ripple 0.1.\n', freqStr, n)
