%##########################################################################
%--------------------------------------------------------------------------
% This function takes in a structure of model parameters as input, solves
% for the stationary equilibrium and returns the outcome.
% The input structure contains the following variables
% N: number of workers
% ce: entry cost denominated in units of labor
% cf: operating cost denominated in units of labor
% ns: length of productivity vector
% svec: productivity vector
% F: transition matrix
% G: startup productivity distribution
% tol: convergence criterion
%-----------------------------
% The output is a structure with the following variables
% pstar: price
% sstar_ind: exit threshold
% mstar: mass of (unmeasured) entrants
% mustar: firm-size distribution before exit stage
% Gnew: startup distribution modified to account for exit
% Fnew: transition matrix modified to account for exit
%--------------------------------------------------------------------------
%##########################################################################

function f = stationary(params)

% Get parameters
cf = params.cf ;
ns = params.ns ;
svec = params.svec ;
F = params.F ;
G = params.G ;
pstar = params.pstar ;
wstar = params.wstar ;
mstar = params.mstar ;

% Solve for value function v and threshold sstar_ind
v = vfi(F, pstar, wstar, cf, params)  ;
sstar_ind = find(v) ;
if (isempty(sstar_ind))
    results.flag = 0 ;
    f = results ;
    return
else
    sstar_ind = sstar_ind(1) ;
end

ce = sum(v.*G);
params.ce = ce ;

% Find stationary distribution

% Start with stationary distribution before exit
Gnew = G ;
Gnew(1:sstar_ind-1) = 0 ;

% New transition matrix to account for exit
Fnew = F ;
Fnew(:,1:sstar_ind-1) = zeros(ns, sstar_ind-1) ;

% Find m0 such that labor markets clear
I = eye(ns) ;
mustar = (I - Fnew) \ (mstar*Gnew) ;
mustar(1:sstar_ind-1) = 0 ;

% calculate employment at each firm
nstar = n(params, svec, pstar, wstar);

nstar_total = nstar + cf;
nstar_total(1:sstar_ind-1) = 0 ;

N = sum(nstar_total.*mustar) + mstar*ce;

avg_stsize_n = sum(nstar_total.*Gnew/sum(Gnew)) ;
avg_fsize_n = sum(nstar_total.*mustar)/sum(mustar) ;
startup_rate = mstar*sum(Gnew)/(sum(mustar)) ;

% Decompose employment
N_prod  = sum(nstar .* mustar) ;
N_fixed = cf * sum(mustar) ;
N_entry = mstar * ce ;

% Output
theta = params.theta ;
Y = sum(pstar * exp(svec) .* nstar.^theta .* mustar) ;

% Average productivity
avg_productivity = Y / N_prod ;

% Exit rate: fraction of incumbents that exit next period
% Firms at s >= sstar_ind can transition to s' < sstar_ind
num_firms = sum(mustar) ;
exit_mass = 0 ;
for i = sstar_ind:ns
    for ip = 1:sstar_ind-1
        exit_mass = exit_mass + F(ip, i) * mustar(i) ;
    end
end
exit_rate = exit_mass / num_firms ;

% Job creation and destruction (continuing firms)
job_creation = 0 ;
job_destruction = 0 ;
for i = sstar_ind:ns
    if mustar(i) > 0
        for ip = 1:ns
            if ip >= sstar_ind
                dn = nstar_total(ip) - nstar_total(i) ;
                job_creation    = job_creation    + max(0,  dn) * F(ip, i) * mustar(i) ;
                job_destruction = job_destruction + max(0, -dn) * F(ip, i) * mustar(i) ;
            else
                % exiting firm destroys all jobs
                job_destruction = job_destruction + nstar_total(i) * F(ip, i) * mustar(i) ;
            end
        end
    end
end

% Entry-related creation
for i = sstar_ind:ns
    if Gnew(i) > 0
        job_creation = job_creation + mstar * Gnew(i) * nstar_total(i) ;
    end
end

job_turnover = (job_creation + job_destruction) / (2 * N_prod) ;

% Store results
results.v = v ;
results.mstar = mstar ;
results.mustar = mustar ;
results.nstar = nstar ;
results.nstar_total = nstar_total ;
results.Gnew = Gnew ;
results.Fnew = Fnew ;
results.sstar_ind = sstar_ind ;
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

f = results ;
