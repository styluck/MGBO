clear; clc; close all;

% % prepare dataset
seed = 40;
randn('state',seed);
a = 5;
r =  a;
n = 2^a;

X = abs(randn(n,r));


% params for BPGM algos
opts_BPGM.gamma = 0.2;
opts_BPGM.rho = 1e-2;

opts_BPGM.maxItr   = 1;
opts_BPGM.record   = 1;
opts_BPGM.mxitr    = 1000;
opts_BPGM.sub_mxitr= 5;
opts_BPGM.ftol     = 1.0e-5;
opts_BPGM.gtol     = 1.0e-10;

opts_BPGM.n        = n;
opts_BPGM.beta     = 1.0;   
opts_BPGM.mu       = 5; %penalty parameter

% main
Zinit=randn(n,r); % initial point
opts_BPGM.r = r;

a = tic;
Bk = orth(Zinit - repmat(sum(Zinit),n,1)/n);

obj_fun = @(X) obj_func_0(X);

% Z = abs(randn(n,r));
% Lamb = (sum(Z).^(-1))';
% obj_fun = @(X) obj_func_svd(X, Z, Lamb);

H = MGBO_BB(Bk, obj_fun, @moreau_hc, @pen_hc, opts_BPGM);

H = sign(H);

cputime = toc(a);

%% calculate feasibility
Obj = obj_fun(H);
feaSt = norm(H'*H - n*eye(r),'fro');
feaKer = norm(H'*ones(n,1), 'fro');
H
fprintf('\nOrth: %2.2f, Ker: %2.2f, Obj:%2.2f\n', feaSt,feaKer, Obj);

% [EOF]
