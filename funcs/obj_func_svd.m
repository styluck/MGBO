function [fval, grad] = obj_func_svd(X, Z, lmb)

ZX = Z'*X;
lZX = lmb.*ZX;
% diagAX = diagA.*X;
fval = -.5*sum(sum(ZX.*lZX));

if nargout > 1
    grad = -( Z*lZX);
end

% fval = .5*(sum(sum(X.*X)) - sum(sum(ZX.*lZX)));
% % fval = .5*(sum(sum(X.*diagAX)) - sum(sum(ZX.*lZX)));
% % 
% if nargout > 1
%     grad = (X - Z*lZX);
% %     grad = (diagAX - Z*lZX);
% end
% [EOF]
