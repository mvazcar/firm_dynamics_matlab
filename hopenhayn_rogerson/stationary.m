%##########################################################################
%--------------------------------------------------------------------------
% This function takes in a structure of model parameters as input, solves
% for the stationary equilibrium and returns the outcome.
% Hopenhayn-Rogerson (1993): 2D distribution over (s, n).
%
% Index convention (matching HR93 notation):
%   i_s      = index for current productivity s
%   i_sprime = index for next-period productivity s'
%   i_n      = index for current employment n
%   i_nprime = index for next-period employment n' (the choice variable)
%   i_nprime2 = second use of n' index (next period's choice, from next period's perspective)
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

% Exit cost for each i_nprime: nn x 1
phi_exit = tau * wstar * nvec ;

% Entry cost: entrants have n = 0 = nvec(1), so phi_exit(1) = 0
% Entrants with i_s stay if v(i_s, 1) >= 0
Gnew = G ;
Gnew(v(:, 1) < 0) = 0 ;

% Free entry: ce = expected value of entering
ce = sum(max(v(:, 1), 0) .* G) ;

if ce <= 0
    results.flag = 0 ;
    f = results ;
    return
end

% Exit flag: firm at (i_sprime, i_nprime) exits if v(i_sprime, i_nprime) < -phi_exit(i_nprime)
exit_flag = v < -phi_exit' ;   % ns x nn, logical

% --- Distribution iteration ---
mu = zeros(ns, nn) ;

for iter = 1:5000
    mu_new = zeros(ns, nn) ;

    % --- Incumbents ---
    for i_n = 1:nn
        for i_s = 1:ns
            if mu(i_s, i_n) > 0
                i_nprime = npol_ind(i_s, i_n) ;          % this period's labor choice
                for i_sprime = 1:ns                       % next period's shock
                    if ~exit_flag(i_sprime, i_nprime)     % stays
                        i_nprime2 = npol_ind(i_sprime, i_nprime) ; % n' from next period's perspective
                        mu_new(i_sprime, i_nprime2) = mu_new(i_sprime, i_nprime2) ...
                            + F(i_sprime, i_s) * mu(i_s, i_n) ;
                    end
                end
            end
        end
    end

    % --- Entrants ---
    % mstar entrants draw s from Gnew, have n_{-1} = 0 = nvec(1)
    for i_s = 1:ns
        if Gnew(i_s) > 0 && ~exit_flag(i_s, 1)
            i_nprime = npol_ind(i_s, 1) ;
            mu_new(i_s, i_nprime) = mu_new(i_s, i_nprime) + mstar * Gnew(i_s) ;
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
    % Entrant at i_s with n=0 chooses nvec(npol_ind(i_s, 1))
    entrant_n = nvec(npol_ind(:, 1)) ;
    avg_stsize_n = sum((entrant_n + cf) .* Gnew) / Gnew_mass ;
else
    avg_stsize_n = 0 ;
end

% Startup rate
startup_rate = mstar * Gnew_mass / num_firms ;

% Exit rate: fraction of firms that exit
% A firm at (i_s, i_n) with choice i_nprime exits next period at (i_sprime, i_nprime)
% with prob F(i_sprime, i_s) when exit_flag(i_sprime, i_nprime) = true
exit_mass = 0 ;
for i_n = 1:nn
    for i_s = 1:ns
        if mustar(i_s, i_n) > 0
            i_nprime = npol_ind(i_s, i_n) ;
            for i_sprime = 1:ns
                if exit_flag(i_sprime, i_nprime)
                    exit_mass = exit_mass + F(i_sprime, i_s) * mustar(i_s, i_n) ;
                end
            end
        end
    end
end
exit_rate = exit_mass / num_firms ;

% Job creation and destruction (cross-period reallocation)
% A firm at (i_s, i_n) chose i_nprime = npol_ind(i_s, i_n) this period.
% Next period it draws i_sprime and, if it survives, chooses npol_ind(i_sprime, i_nprime).
% Turnover = |nvec(npol_ind(i_sprime, i_nprime)) - nvec(i_nprime)| weighted by transition prob and mu.
job_creation = 0 ;
job_destruction = 0 ;
for i_n = 1:nn
    for i_s = 1:ns
        if mustar(i_s, i_n) > 0
            i_nprime = npol_ind(i_s, i_n) ;            % this period's labor choice
            n_now = nvec(i_nprime) ;
            for i_sprime = 1:ns
                if ~exit_flag(i_sprime, i_nprime)       % firm survives
                    i_nprime2 = npol_ind(i_sprime, i_nprime) ; % n' from next period's perspective
                    dn = nvec(i_nprime2) - n_now ;
                    job_creation    = job_creation    + max(0,  dn) * F(i_sprime, i_s) * mustar(i_s, i_n) ;
                    job_destruction = job_destruction + max(0, -dn) * F(i_sprime, i_s) * mustar(i_s, i_n) ;
                else
                    % exiting firm destroys all workers
                    job_destruction = job_destruction + n_now * F(i_sprime, i_s) * mustar(i_s, i_n) ;
                end
            end
        end
    end
end

% Add entry-related creation
for i_s = 1:ns
    if Gnew(i_s) > 0 && ~exit_flag(i_s, 1)
        i_nprime = npol_ind(i_s, 1) ;
        job_creation = job_creation + mstar * Gnew(i_s) * nvec(i_nprime) ;
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
