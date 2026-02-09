% This function takes model parameters as inputs, does the value function
% iteration over (s, n) and returns the value function and policy indices.

function [f, npol_ind] = vfi(F, pstar, wstar, cf, tau, params)
beta = params.beta ;
svec = params.svec ;
nvec = params.nvec ;
ns = params.ns ;
nn = params.nn ;
tol = params.tol ;
theta = params.theta ;

% Precompute flow profit for each (s_i, nprime_k) pair: ns x nn
pi_mat = pstar * exp(svec) .* nvec'.^theta - wstar * nvec' - wstar * cf ;

% Firing cost for each (nprime_k, n_j) pair: nn x nn
% phi_mat(k,j) = cost of going FROM nvec(j) TO nvec(k)
phi_mat = tau * wstar * max(0, nvec - nvec') ;   % nn x nn

% Exit cost for each nprime_k: nn x 1
% phi_exit(k) = cost of exiting next period with nvec(k) workers
phi_exit = tau * wstar * nvec ;                   % nn x 1

% Initialize
v = zeros(ns, nn) ;
Tv = ones(ns, nn) ;
npol_ind = ones(ns, nn) ;
maxiter = 1000 ;
iter = 0 ;

while max(abs(Tv(:) - v(:))) > tol && iter <= maxiter
    v = Tv ;

    % Next-period value with exit option: max(v(s', k), -phi_exit(k))
    vmax = max(v, -phi_exit') ;                    % ns x nn

    % Expected continuation: F' * vmax
    Ev = F' * vmax ;                                % ns x nn

    % Maximize over nprime choice k, for each state (s_i, n_j)
    for j = 1:nn
        % obj(i, k) = pi_mat(i,k) - phi_mat(k,j) + beta * Ev(i,k)
        obj = pi_mat - phi_mat(:, j)' + beta * Ev ;    % ns x nn
        [Tv(:, j), npol_ind(:, j)] = max(obj, [], 2) ; % max over k
    end

    iter = iter + 1 ;
end

f = v ;
