function [fval, grad] = obj_func_0(X)


% diagAX = diagA.*X;
fval = 0;

if nargout > 1
    grad = 0;
end

% fval = .5*(sum(sum(X.*X)) - sum(sum(ZX.*lZX)));
% % fval = .5*(sum(sum(X.*diagAX)) - sum(sum(ZX.*lZX)));
% % 
% if nargout > 1
%     grad = (X - Z*lZX);
% %     grad = (diagAX - Z*lZX);
% end
% [EOF]