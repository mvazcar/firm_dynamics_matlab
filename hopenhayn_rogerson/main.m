% =========================================================================
% Firm Dynamics Model - Labor Only
% Main: Hopenhayn-Rogerson (1993)
% =========================================================================
%% -------------------------------------------------------------------------
% Initialize parameters and other model primitives
% -------------------------------------------------------------------------
clc
clear
close all

beta  = 1/1.04 ; % Discount rate
theta   = 0.64   ; % Returns-to-scale
ns    = 100    ; % Number of points in productivity grid
tol   = 1e-8   ; % Specify tolerance

pstar = 1 ;   % Equilibrium output price
wstar = 1;
mstar = 100 ; % Mass of entrants in steady state

% Create a structure that we want to transfer to functions
params.beta = beta ;
params.theta = theta ;
params.tol = tol ;
params.ns = ns ;
params.pstar = pstar ;
params.wstar = wstar ;
params.mstar = mstar ;

% Set parameters
rho       = 0.984150757243253;    % Persistence  of AR1
mu        = -1.436111629482697;   % Long-run AR1 mean
sigma     = 0.245520815536363;    % Standard deviation of AR1 innovations
cf        = 24.308026243791222;   % Fixed cost

mu0       = -4.344376541584754;   % Mean of startup distribution (\mu_g)
sigma0    = 1.331137767741511;    % Standard deviation of startup distribution (\sigma_G^2)

%% ------------------------------------------------------------------------

params.cf = cf ;
params.sigma0 = sigma0 ;

[svec, F] = tauchen(mu, rho, sigma, ns) ;
F = transpose(F) ; % Transpose so that columns sum to one

% Add to params
params.svec = svec ;
params.F = F ;

% Startup productivity distribution
step = svec(2) - svec(1) ;
G = normcdf(svec + step/2, mu0, sigma0) ;
G(2:end-1) = G(2:end-1) - G(1:end-2) ;
G(end) = 1 - sum(G(1:end-1)) ;
params.G = G ;

%% -------------------------------------------------------------------------
% Employment grid
% -------------------------------------------------------------------------
nn = 300 ;
n_max = 2 * max(n(params, svec, pstar, wstar)) ;
nvec = [0; logspace(log10(1), log10(n_max), nn-1)'] ;

params.nn = nn ;
params.nvec = nvec ;

%% -------------------------------------------------------------------------
% Solve for stationary equilibrium at different tau values
% -------------------------------------------------------------------------

% Benchmark: tau = 0
params.tau = 0 ;
fprintf('Solving tau = 0 ...\n') ;
results_0 = stationary(params) ;

% tau = 0.1
params.tau = 0.1 ;
fprintf('Solving tau = 0.1 ...\n') ;
results_1 = stationary(params) ;

% tau = 0.5
params.tau = 0.5 ;
fprintf('Solving tau = 0.5 ...\n') ;
results_2 = stationary(params) ;

%% -------------------------------------------------------------------------
% Display results (Table 3 style)
% -------------------------------------------------------------------------
fprintf('\n===== Hopenhayn-Rogerson (1993) Comparative Statics =====\n') ;
fprintf('%-25s %12s %12s %12s\n', '', 'tau=0', 'tau=0.1', 'tau=0.5') ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Entry cost (ce)', results_0.ce, results_1.ce, results_2.ce) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Exit rate', results_0.exit_rate, results_1.exit_rate, results_2.exit_rate) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Startup rate', results_0.startup_rate, results_1.startup_rate, results_2.startup_rate) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Avg firm size', results_0.avg_fsize_n, results_1.avg_fsize_n, results_2.avg_fsize_n) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Avg startup size', results_0.avg_stsize_n, results_1.avg_stsize_n, results_2.avg_stsize_n) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Avg productivity', results_0.avg_productivity, results_1.avg_productivity, results_2.avg_productivity) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Job turnover', results_0.job_turnover, results_1.job_turnover, results_2.job_turnover) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Total employment', results_0.N, results_1.N, results_2.N) ;
fprintf('%-25s %12.4f %12.4f %12.4f\n', 'Output', results_0.Y, results_1.Y, results_2.Y) ;
