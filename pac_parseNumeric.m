function numbers = pac_parseNumeric(inputText)
%PAC_PARSENUMERIC Parse a numeric GUI field without eval or str2num.
% Allowed syntax: decimal, scientific notation and whitespace/comma pairs.
% Empty input returns [] for an optional argument.
%
% Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.

if isempty(inputText)
    numbers=[];
    return
end
if ~(ischar(inputText) || (isstring(inputText) && isscalar(inputText)))
    error('PACT numeric input must be text.');
end
s=strtrim(char(inputText));
if isempty(s)
    numbers=[];
    return
end
s=strrep(s,'[',' ');
s=strrep(s,']',' ');
s=strrep(s,',',' ');
parts=regexp(strtrim(s),'\s+','split');
numbers=zeros(1,numel(parts));
for k=1:numel(parts)
    if isempty(regexp(parts{k},'^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$','once'))
        error('PACT numeric inputs must be literal numbers, not expressions.');
    end
    numbers(k)=str2double(parts{k});
end
if any(~isfinite(numbers))
    error('PACT numeric input must be finite.');
end
end
