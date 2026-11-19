function [order, df] = pac_firwsOrder(srate, cutoffs, requestedDf)
%PAC_FIRWSORDER Hamming windowed-sinc order whose transition fits the band.
%   [order, df] = pac_firwsOrder(srate, [low high], requestedDf)
% Returns the even firfilt/pop_firws filter order for a Hamming windowed-sinc
% band-pass filter with -6 dB cutoffs at [low high]. The transition width is
% requestedDf, capped at pop_firws' limit min(2*low, 2*(Nyquist-high),
% high-low). Above that limit the two transition bands overlap, the centre
% gain drops below 1 and pop_firws warns "Filter order too low".
% df is the transition width actually achieved by the even order.
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

if ~isscalar(srate) || srate <= 0 || numel(cutoffs) ~= 2 || cutoffs(1) <= 0 || ...
        cutoffs(2) <= cutoffs(1) || cutoffs(2) >= srate/2 || ~isscalar(requestedDf) || requestedDf <= 0
    error('PACT:InvalidInput','Invalid sampling rate, cutoffs or transition width.');
end
nyquist = srate/2;
maxDf = min([2*cutoffs(1), 2*cutoffs(2), 2*(nyquist-cutoffs(1)), 2*(nyquist-cutoffs(2)), diff(cutoffs)]);
target = min(requestedDf, maxDf);
order = ceil(3.3/(target/srate)/2)*2;   % Hamming: df = 3.3*srate/order
df = 3.3*srate/order;
end
