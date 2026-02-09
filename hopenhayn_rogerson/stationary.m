%##########################################################################
%--------------------------------------------------------------------------
% This function takes in a structure of model parameters as input, solves
% for the stationary equilibrium and returns the outcome.
% Hopenhayn-Rogerson (1993): 2D distribution over (s, n).
%--------------------------------------------------------------------------
%##########################################################################

function f = stationary(params)

% Get parameters
cf = params.cf ;
ns = params.ns ;
nn = params.nn ;
svec = params.svec ;
nvec = params.nvec ;
F = params.F ;
G = params.G ;
pstar = params.pstar ;
wstar = params.wstar ;
mstar = params.mstar ;
tau = params.tau ;
theta = params.theta ;

% Solve for value function v and policy npol_ind
[v, npol_ind] = vfi(F, pstar, wstar, cf, tau, params) ;

% Exit cost for each nprime_k: nn x 1
phi_exit = tau * wstar * nvec ;

% Entry cost: entrants have n = 0 = nvec(1), so phi_exit(1) = 0
% Entrants with s_i stay if v(i, 1) >= 0
Gnew = G ;
Gnew(v(:, 1) < 0) = 0 ;

% Free entry: ce = expected value of entering
ce = sum(max(v(:, 1), 0) .* G) ;

if ce <= 0
    results.flag = 0 ;
    f = results ;
    return
end

% Exit flag: firm at (s', nprime_k) exits if v(s', k) < -phi_exit(k)
exit_flag = v < -phi_exit' ;   % ns x nn, logical

% --- Distribution iteration ---
mu = zeros(ns, nn) ;

for iter = 1:5000
    mu_new = zeros(ns, nn) ;

    % --- Incumbents ---
    for j = 1:nn
        for i = 1:ns
            if mu(i, j) > 0
                k = npol_ind(i, j) ;          % this period's choice
                for ip = 1:ns                  % next period's shock
                    if ~exit_flag(ip, k)       % stays
                        kp = npol_ind(ip, k) ; % next period's choice
                        mu_new(ip, kp) = mu_new(ip, kp) + F(ip, i) * mu(i, j) ;
                    end
                end
            end
        end
    end

    % --- Entrants ---
    % mstar entrants draw s from Gnew, have n_{-1} = 0 = nvec(1)
    for i = 1:ns
        if Gnew(i) > 0 && ~exit_flag(i, 1)
            k = npol_ind(i, 1) ;
            mu_new(i, k) = mu_new(i, k) + mstar * Gnew(i) ;
        end
    end

    if max(abs(mu_new(:) - mu(:))) < 1e-6
        break
    end
    mu = mu_new ;
end

mustar = mu ;

% --- Aggregates ---

% Employment
N_prod = sum(nvec' .* mustar, 'all') ;          % production workers
N_fixed = cf * sum(mustar, 'all') ;               % fixed cost workers
N_entry = mstar * ce ;                             % entry labor
N = N_prod + N_fixed + N_entry ;

% Output
Y = sum((pstar * exp(svec) .* nvec'.^theta) .* mustar, 'all') ;

% Average productivity
avg_productivity = Y / N_prod ;

% Average firm size (production workers)
num_firms = sum(mustar, 'all') ;
avg_fsize_n = (N_prod + cf * num_firms) / num_firms ;

% Average startup size
Gnew_mass = sum(Gnew) ;
if Gnew_mass > 0
    % Entrant at s_i with n=0 chooses nvec(npol_ind(i,1))
    entrant_n = nvec(npol_ind(:, 1)) ;
    avg_stsize_n = sum((entrant_n + cf) .* Gnew) / Gnew_mass ;
else
    avg_stsize_n = 0 ;
end

% Startup rate
startup_rate = mstar * Gnew_mass / num_firms ;

% Exit rate: fraction of firms that exit
% A firm at (s_i, n_j) exits next period at (s', k) with prob sum of F(ip,i)*exit_flag(ip,k)
exit_mass = 0 ;
for j = 1:nn
    for i = 1:ns
        if mustar(i, j) > 0
            k = npol_ind(i, j) ;
            for ip = 1:ns
                if exit_flag(ip, k)
                    exit_mass = exit_mass + F(ip, i) * mustar(i, j) ;
                end
            end
        end
    end
end
exit_rate = exit_mass / num_firms ;

% Job creation and destruction (cross-period reallocation)
% A firm at (s_i, n_j) chose nprime = nvec(k) this period.
% Next period it draws s' and, if it survives, chooses nvec(npol_ind(s', k)).
% Turnover = |nvec(npol_ind(s', k)) - nvec(k)| weighted by F(s'|s) and mu.
job_creation = 0 ;
job_destruction = 0 ;
for j = 1:nn
    for i = 1:ns
        if mustar(i, j) > 0
            k = npol_ind(i, j) ;           % this period's employment index
            n_now = nvec(k) ;
            for ip = 1:ns
                if ~exit_flag(ip, k)        % firm survives
                    kp = npol_ind(ip, k) ;  % next period's employment index
                    dn = nvec(kp) - n_now ;
                    job_creation    = job_creation    + max(0,  dn) * F(ip, i) * mustar(i, j) ;
                    job_destruction = job_destruction + max(0, -dn) * F(ip, i) * mustar(i, j) ;
                else
                    % exiting firm destroys all workers
                    job_destruction = job_destruction + n_now * F(ip, i) * mustar(i, j) ;
                end
            end
        end
    end
end

% Add entry-related creation
for i = 1:ns
    if Gnew(i) > 0 && ~exit_flag(i, 1)
        k = npol_ind(i, 1) ;
        job_creation = job_creation + mstar * Gnew(i) * nvec(k) ;
    end
end

job_turnover = (job_creation + job_destruction) / (2 * N_prod) ;

% Store results
results.v = v ;
results.npol_ind = npol_ind ;
results.mstar = mstar ;
results.mustar = mustar ;
results.Gnew = Gnew ;
results.ce = ce ;
results.N = N ;
results.N_prod = N_prod ;
results.N_fixed = N_fixed ;
results.N_entry = N_entry ;
results.Y = Y ;
results.avg_productivity = avg_productivity ;
results.avg_fsize_n = avg_fsize_n ;
results.avg_stsize_n = avg_stsize_n ;
results.startup_rate = startup_rate ;
results.exit_rate = exit_rate ;
results.job_creation = job_creation ;
results.job_destruction = job_destruction ;
results.job_turnover = job_turnover ;
results.tau = tau ;

f = results ;
