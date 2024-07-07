function [f, g] = moreau_hc(X, sqn, mu, gmma)
 
% global sqn mu gmma;

gmma_mu = mu/gmma; % 1/beta*mu

sgnX = sign(X);

% calculate Z = prox(X)
Z = sign(X)*sqn;
absX = abs(X);
set1 = absX < sqn;
set2 = absX > sqn + gmma;

Z(set1) = X(set1);
Z(set2) = X(set2) - gmma*sign(X(set2));

XZ = X - Z;

% f = env(X- 1/beta*Lambda)
f = mu*sum(sum(max(0, -Z - sqn) + max(0, Z - sqn))) ... % \rho h(Z)
    + (0.5*gmma_mu)*norm(XZ, 'fro')^2; % + \rho/\gamma \|X-Z\|^2

if nargout > 1
     g = gmma_mu*XZ;
end

end
% [EOF]

