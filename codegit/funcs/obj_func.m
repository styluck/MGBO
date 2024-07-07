function [fval, grad] = obj_func(X, H, varargin)

HX = H*X;

fval = .5*sum(sum(X.*HX));

if nargout > 1
    grad = HX;
end

% [EOF]