%--------------------------------------------------------------------------
% Parse varargin Name-Value pairs
%--------------------------------------------------------------------------

function [found, Value, vararg] = argParse(vararg, Name)

NameInd = cellfun(@(x) isequal(lower(string(Name)),x),cellfun(@lower,vararg,'un',0));

% NameInd = cellfun(@(x) (ischar(x) || (isstring(x) && length(x)==1)) && strcmpi(x, Name), vararg);

if nnz(NameInd) == 0
    found = false;
    Value = [];
elseif nnz(NameInd) == 1
    found = true;
    jj = find(NameInd);
    Value = vararg{jj+1};
    vararg([jj jj+1]) = [];
else
    error('found multiple instances of %s argument',upper(Name));
end
end