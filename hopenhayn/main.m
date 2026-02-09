% =========================================================================
% Firm Dynamics Model - Labor Only
% Main: Hopenhayn (1992) with Hopenhayn, Neira, Singhania (2022)
% =========================================================================
%% -------------------------------------------------------------------------
% Initialize parameters and other model primitives
% -------------------------------------------------------------------------
clc
clear
close all

beta  = 1/1.04 ; % Discount rate
theta   = 0.64   ; % Returns-to-scale
ns    = 100    ; % Number of points in grid
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
F = transpose(F) ; % Transpose so that column sum to one

% Add to params
params.svec = svec ;
params.F = F ;

% Startup productivity distribution
step = svec(2) - svec(1) ;
G = normcdf(svec + step/2, mu0, sigma0) ;
G(2:end-1) = G(2:end-1) - G(1:end-2) ;
G(end) = 1 - sum(G(1:end-1)) ;
params.G = G ;

% Solve for stationary equilibrium in 1940, store results in stateqbm
results = stationary(params) ;

%% -------------------------------------------------------------------------
% Display results
% -------------------------------------------------------------------------
fprintf('\n===== Hopenhayn (1992) Stationary Equilibrium =====\n') ;
fprintf('%-25s %12.4f\n', 'Entry cost (ce)', results.ce) ;
fprintf('%-25s %12.4f\n', 'Exit rate', results.exit_rate) ;
fprintf('%-25s %12.4f\n', 'Startup rate', results.startup_rate) ;
fprintf('%-25s %12.4f\n', 'Avg firm size', results.avg_fsize_n) ;
fprintf('%-25s %12.4f\n', 'Avg startup size', results.avg_stsize_n) ;
fprintf('%-25s %12.4f\n', 'Avg productivity', results.avg_productivity) ;
fprintf('%-25s %12.4f\n', 'Job turnover', results.job_turnover) ;
fprintf('%-25s %12.4f\n', 'Total employment', results.N) ;
fprintf('%-25s %12.4f\n', 'Output', results.Y) ;