%%
clear all; clc;
randn('state', 0);
rand('state', 0);

n = 100;

% generate a well-conditioned positive definite matrix
% (for faster convergence)
P = rand(n);
P = P + P';
[V D] = eig(P);
P = V*diag(1+rand(n,1))*V';

q = randn(n,1);
r = 0;

l = randn(n,1);
u = randn(n,1);
lb = min(l,u);
ub = max(l,u);

% ADMM
tstart=tic
[x history] = quad_ADMM(P, q, r, lb, ub, 1.0, 1.0);
toc(tstart)

% normal algirithm: interior point convex
tstart = tic
[x,fval] = quadprog(P, q, [], [], [], [], lb, ub);
toc(tstart)

% linear programming
f = randi([-100,100], 1, n);
tic
[x fval] = linprog(f,[],[],[],[],lb,ub);
toc

%% portfolio optimization problem
clear all; clc;
load('port5.mat','Correlation','stdDev_return','mean_return')
% Calculate covariance matrix from correlation matrix.
Covariance = Correlation .* (stdDev_return * stdDev_return');
nAssets = numel(mean_return); r = 0.002;     % number of assets and desired return
Aeq = ones(1,nAssets); beq = 1;              % equality Aeq*x = beq
Aineq = -mean_return'; bineq = -r;           % inequality Aineq*x <= bineq
lb = zeros(nAssets,1); ub = ones(nAssets,1); % bounds lb <= x <= ub
c = zeros(nAssets,1);                        % objective has no linear term; set it to zero
options = optimoptions('quadprog','Algorithm','interior-point-convex');
% Set additional options: turn on iterative display, and set a tighter optimality termination tolerance.
options = optimoptions(options,'Display','iter','TolFun',1e-10);

% Call solver and measure wall-clock time.
tic
[x1,fval1] = quadprog(Covariance,c,Aineq,bineq,Aeq,beq,lb,ub,[],options);
toc

% Plot results.
%plotPortfDemoStandardModel(x1)

% linear programming
f = randi([-100,100],1,nAssets);
options2 = optimoptions('linprog','Algorithm','interior-point')
tic
[xlinear, fval2] = linprog(f, Aineq, bineq, Aeq, beq, lb, ub, [], options2);
toc