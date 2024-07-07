function h = pen_hc(x, sqn)

if nargin < 2 || isempty(sqn)
    sqn = 1/sqrt(size(x,1));
end

h = sum(sum(max(0,-x-sqn) + max(0,x-sqn)));
end

% [EOF]