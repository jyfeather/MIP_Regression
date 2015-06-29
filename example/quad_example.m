randn('state', 0);
rand('state', 0);

n = 100;

% generate a well-conditioned positive definite matrix
% (for faster convergence)
P = rand(n);
P = P + P';
[V D] = eig(P);
P = V*diag(1+rand(n,1))*V';
P = (P+P')/2

q = randn(n,1);
r = 0;

l = randn(n,1);
u = randn(n,1);
lb = min(l,u);
ub = max(l,u);

%% ADMM
tstart=tic
[x history] = quad_ADMM(P, q, r, lb, ub, 1.0, 1.0);
toc(tstart)

%% normal algirithm: interior point convex
tstart = tic
[x,fval] = quadprog(P, q, [], [], [], [], lb, ub)
toc(tstart)