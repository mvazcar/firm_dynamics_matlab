% This function takes model parameters as inputs, does the value function
% iteration over (s, n) and returns the value function and policy indices.
%
% Index convention (matching HR93 notation):
%   i_s      = index for current productivity s
%   i_sprime = index for next-period productivity s'
%   i_n      = index for current employment n
%   i_nprime = index for next-period employment n' (the choice variable)
%
% Grids:
%   svec(i_s)      = productivity level       (ns x 1)
%   nvec(i_n)      = employment level          (nn x 1)
%
% Matrices returned:
%   v(i_s, i_n)        = value function
%   npol_ind(i_s, i_n) = index of optimal n' given state (s, n)

function [f, npol_ind] = vfi(F, pstar, wstar, cf, tau, params)
beta = params.beta ;
svec = params.svec ;
nvec = params.nvec ;
ns = params.ns ;
nn = params.nn ;
tol = params.tol ;
theta = params.theta ;

% Precompute flow profit for each (i_s, i_nprime) pair: ns x nn
% pi(s, nprime) = p * exp(s) * nprime^theta - w * nprime - w * cf
pi_mat = pstar * exp(svec) .* nvec'.^theta - wstar * nvec' - wstar * cf ;

% Firing cost: phi(n, nprime) = tau * w * max(0, n - nprime) as in HR93
% Broadcasting: nvec (nn x 1) = rows = n, nvec' (1 x nn) = cols = nprime
% phi_mat(i_n, i_nprime) = tau * w * max(0, n - nprime)
phi_mat = tau * wstar * max(0, nvec - nvec') ;   % nn x nn

% Exit cost: firing all workers when exiting
% phi_exit(i_nprime) = tau * w * nvec(i_nprime)
phi_exit = tau * wstar * nvec ;                   % nn x 1

% Initialize
v = zeros(ns, nn) ;
Tv = ones(ns, nn) ;
npol_ind = ones(ns, nn) ;
maxiter = 1000 ;
iter = 0 ;

while max(abs(Tv(:) - v(:))) > tol && iter <= maxiter
    v = Tv ;

    % Exit option: max(v(i_sprime, i_nprime), -phi_exit(i_nprime))
    vmax = max(v, -phi_exit') ;                    % ns x nn

    % Expected continuation: E[vmax | s] = F' * vmax
    Ev = F' * vmax ;                                % ns x nn

    % Maximize over n' for each state (s, n)
    for i_n = 1:nn
        % obj(i_s, i_nprime) = pi(s, nprime) - phi(n, nprime) + beta * E[vmax(s, nprime)]
        obj = pi_mat - repmat(phi_mat(i_n, :), ns, 1) + beta * Ev ;    % ns x nn, note the use of repmat to broadcast
        [Tv(:, i_n), npol_ind(:, i_n)] = max(obj, [], 2) ; % max over i_nprime
    end

    iter = iter + 1 ;
end

% -----------------------------------------------------------------------
% Value function iteration: fill Tv(ns x nn) one column at a time.
%
% For each current employment i_n:
%   1. Build obj(i_s, i_nprime) = profit(s,n') - firing_cost(n,n') + beta*E[v(s,n')]
%      - pi_mat is ns x nn (depends on s and n')
%      - phi_mat(i_n,:) is 1 x nn, broadcast (repeated) ns times to ns x nn
%      - beta*Ev is ns x nn (depends on s and n')
%
%   2. max(obj, [], 2) maximizes over columns (n') for each row (s),
%      returning an ns x 1 vector: the best value for each s given this n.
%
%   3. That ns x 1 vector is stored in column i_n of Tv:
%        Tv(:, i_n) = max(obj, [], 2)
%      and the argmax gives the optimal n' index:
%        npol_ind(:, i_n) = argmax over n'
%
% After looping over i_n = 1, ..., nn, the full ns x nn matrix Tv is filled.
% Each entry Tv(i_s, i_n) = maximized value at state (s, n).
% Each entry npol_ind(i_s, i_n) = index of optimal n' at state (s, n).
%
% Slow loop equivalent (for reference):
% for i_n = 1:nn
%   for i_s = 1:ns
%     for i_nprime = 1:nn
%       obj(i_s, i_nprime) = pi_mat(i_s, i_nprime) - phi_mat(i_n, i_nprime) + beta * Ev(i_s, i_nprime) ;
%     end
%   [Tv(i_s, i_n), npol_ind(i_s, i_n)] = max(obj(i_s, :)) ;
%   end
% end
% -----------------------------------------------------------------------


f = v ;
