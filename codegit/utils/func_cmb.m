
function [f, g] = func_cmb(ffun,hfun,varargin)
% combine f function and h function
 [f1, g1] = ffun(varargin{:});
 [f2, g2] = hfun(varargin{:});
 f = f1 + f2;
 g = g1 + g2;

%  ffun = @(X) objQAP_QUAD(X,A,B);
% hfun = @(X) pen_quad(X,mu,lmb);
% [outputs{:}] = ffun;
% hh = @(X) ffun(X) + hfun(X);

% [EOF]
