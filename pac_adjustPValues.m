function adjusted = pac_adjustPValues(pvals, method)
%PAC_ADJUSTPVALUES Adjust a family of p-values across channels.
%   adjusted = pac_adjustPValues(pvals, 'none'|'bonferroni'|'holm'|'fdr')
%   Missing (NaN) p-values are excluded from the testing family.
%   'fdr' implements exact Benjamini-Hochberg adjusted p-values.
%   All returned finite values are in [0,1].
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

if nargin ~= 2
    error('pac_adjustPValues requires pvals and method.');
end
adjusted = nan(size(pvals));
valid = isfinite(pvals);
if any(pvals(valid) < 0 | pvals(valid) > 1)
    error('P-values must lie between 0 and 1.');
end
indices = find(valid);
m = numel(indices);
if m == 0
    return
end
[sorted, order] = sort(pvals(indices));
sorted = sorted(:);
switch lower(method)
    case 'none'
        corrected = sorted;
    case 'bonferroni'
        corrected = min(1, m .* sorted);
    case 'holm'
        corrected = min(1, cummax(sorted .* (m:-1:1)'));
    case 'fdr'
        corrected = min(1, flipud(cummin(flipud(sorted .* m ./ (1:m)'))));
    otherwise
        error('Unknown multiple-comparison method: %s', method);
end
adjusted(indices(order)) = corrected;
end
