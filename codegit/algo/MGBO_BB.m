function [Xk, outs] = MGBO_BB(Xk, ffun, gfun, penfun, opts, varargin)
%[Xnew,iter,objf] = MGBO_BB(Xk,retr,pars,OPTIONS,A,B,lambda)

if ~isfield(opts, 'gtol');      opts.gtol = 1e-3; end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-6; end
if ~isfield(opts, 'xtol');      opts.xtol = 1e-6; end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'tau');      opts.tau  = 1.0e-2; end % stepsize
if ~isfield(opts, 'mintau');   opts.mintau  = 1e-12; end % mini-stepsize
if ~isfield(opts, 'maxtau');   opts.maxtau  = opts.tau; end % max-stepsize
if ~isfield(opts, 'rhols');    opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');      opts.eta  = 0.1; end
if ~isfield(opts, 'gamma');    opts.gamma  = 0.85; end
if ~isfield(opts, 'mxitr');    opts.mxitr  = 100; end
if ~isfield(opts, 'record');   opts.record = 0; end
if ~isfield(opts, 'retr');     opts.retr = 3; end
if ~isfield(opts, 'mu');       opts.mu   = 1; end
if ~isfield(opts, 'manifold'); opts.manifold = 1; end

%%
%-------------------------------------------------------------------------------
% copy parameters
n      = opts.n;
r      = opts.r;
sqtn   = sqrt(n);
gtol   = opts.gtol*sqtn; %
ftol   = opts.ftol*sqtn; %
xtol   = opts.xtol*sqtn; %
rhols  = opts.rhols;
eta    = opts.eta;
tau    = opts.tau; % initial step size
mintau = opts.mintau;
maxtau = opts.maxtau;
record = opts.record; %
kappa  = 5;

mu     = opts.mu;
gmma   = opts.gamma;

global sqn;
sqn = 1/sqtn;

% function initialization
ffun = @(X) ffun(X, varargin{:});
gfunc = @(X) gfun(X, sqn, mu, gmma);

func = @(X) func_cmb(ffun, gfunc, X);

distfunc = @(X) max(max(abs(X - sign(X)*sqn)));

if opts.manifold == 1 %\bnr
    projT = @(X, V) projTB(X, V, n);
elseif opts.manifold == 2 %\snr
    projT = @(X, V) projM(X, V);
end
%% ****************************************************************
% [f1, g1] = ffun(Xk,varargin{:});
% [f2, g2] = gfun(Xk,sqn, mu, gmma);
% Fold = f1 + f2;
% grad = g1 + g2;
[Fold, grad] = func(Xk);
grad = projT(Xk, grad);

Fk_list = zeros(opts.mxitr,1);

outs.nfe = 0;

if (record)
    fprintf('\n itr   nls      g(x)   \t\t  fval  \t   ||x-proj(x)||    normVk    time');
    tstart = tic;
end

Cval = Fold;

for itr = 1:opts.mxitr
    
    nls = 1;
    
    while 1   %% to search a desired step-size
        
        Vk = -tau*grad;
        
        normVk = norm(Vk,'fro');
        sqnormVk = normVk^2;
        
        Xnew = myQRretr(Xk,Vk);
        
        %         [f1, g1] = ffun(Xnew, varargin{:});
        %         [f2, g2] = gfun(Xnew, sqn, mu, gmma);
        %         Fnew = f1 + f2;
        %         grad_new = g1 + g2;
        [Fnew, grad_new] = func(Xnew);
        if (Fnew<= Cval-(0.5/tau)*rhols*sqnormVk)||(nls>=5)
            break;
        end
        
        tau = eta*tau;
        
        nls = nls + 1;
    end
    outs.nfe = outs.nfe + nls;
    
    Fk_list(itr) = Fnew;
    
    % fprintf('nabla f:%2.2f', norm(grad_new,'fro'));
    grad_new = projT(Xnew, grad_new);
    
    norm_grad = norm(grad_new,'fro');
    dist1 = distfunc(Xnew);
    if norm_grad <=gtol && Fnew<=Fold && dist1 < xtol
        break;
        
    elseif (Fnew<=Fold && itr>=5 ...
            && abs(Fnew-Fold)/(1+abs(Fnew))<=ftol) && dist1 < xtol
        break;
    end
    
    %% *************** to estimate the step-size via BB *****************
    
    DetaX = Xnew - Xk;  DetaY = grad_new - grad;
    
    DetaXY = abs(sum(dot(DetaX,DetaY)));
    
    tau1 = norm(DetaX,'fro')^2/DetaXY;
    
    tau2 = DetaXY/norm(DetaY,'fro')^2;
    
    tau = max(min(min(tau1, tau2),maxtau), mintau);
    
    % recorder
    if (record) && rem(itr,10) == 0
        ttime = toc(tstart);
        gapk = penfun(Xk, sqn);
        fk = full(ffun(Xk));
        fprintf('\n %4d  %3d     %3.2e     %5.4e    %6.2e    %6.2e    %2.1f',...
            itr, nls, gapk, fk, dist1, norm_grad, ttime);
    end
        
    Xk = Xnew;  grad = grad_new;
    
    if itr <= kappa
        Cval = max(Fk_list(1:itr));
    else
        Cval = max(Fk_list(itr-kappa+1:itr));
    end
end

% recorder
if (record)
    ttime = toc(tstart);
    gapk = penfun(Xk, sqn);
    fk = full(ffun(Xk));
    fprintf('\n %4d  %3d     %3.2e     %5.4e    %6.2e    %6.2e    %2.1f',...
        itr, nls, gapk, fk, dist1, norm_grad, ttime);
else
    gapk = penfun(Xk, sqn);
end
outs.itr = itr;
outs.fval = Fnew;
outs.gapk = gapk;
%-------------------------------------------------------------------------------
    function Q = myQRretr(Xk, Vk)
        % retraction onto the manifold based on QR decomposition
        mtd = opts.retr;
        switch mtd
            case 1 % exponential map
                XkVk = Xk'*Vk;
                [QQ,R]=qr(-Vk+Xk*XkVk,0);
                Z = zeros(r);
                Q = [Xk, QQ] * expm([ -XkVk, -R'; R, Z])*[eye(r); Z];
                
            case 2 % polar decomposition
                %                 fprintf('polar\n');
                XX = Xk + Vk;
                [u, ~, v] = svd(XX, 'econ'); % #ok
                Q = u*v';
                
            case 3 % qr decomposition
                %                 fprintf('qr\n');
                [Q, RR] = qr(Xk + Vk, 0);
                diagRR = sign(diag(RR));
                ndr = diagRR < 0;
                if nnz(ndr)> 0
                    Q = Q*spdiags(diagRR,0, r, r);
                end
            case 4 % Cayley
                %                 fprintf('Cay\n');
                U =  [Vk, Xk];    V = [Xk, -Vk];       VU = V'*U;
                
                VX = V'*Xk;
                [aa, ~] = linsolve(eye(2*r) + (0.5)*VU, VX);
                Q = Xk - U*(aa);
        end
    end


end

%-------------------------------------------------------------------------------
% projection onto the tangent space of manifold
function U = projTB(Xk, Z, n)
% projection onto the tangent space of the B manifold
Ztilde = Z - repmat(sum(Z),n,1)/n;
XZ = Xk'*Ztilde;

U = Ztilde - .5*Xk*(XZ+XZ');
end

function U = projM(Xk, Ztilde)
% projection onto the tangent space of the Stiefel manifold
XZ = Xk'*Ztilde;

U = Ztilde - .5*Xk*(XZ+XZ');
end

% [EOF]
