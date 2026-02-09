function f = phi(np, n, tau, w)
% Firing cost: tau per job destroyed, in wage units
f = tau * w * max(0, n - np) ;
