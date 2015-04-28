function [ out ] = mbe_iff(cond, val_true, val_false)
%IFF Ternary conditinal operator, like "cond ? a:b" in C language. 
%   Usage: iff(condition, value_if_true, value_if_false)
% Though, not 100% efficient as in C, since here we force evaluating both sides of
% the condition.
    if (cond)
        out = val_true;
    else
        out = val_false;
    end
end

