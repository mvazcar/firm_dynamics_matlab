% ##########################################################################
% Optimal Labor from Firm FOC - Static
% ##########################################################################

function f = n(params, s, p, w)
theta = params.theta ;
f = ((p*exp(s)*theta)/w).^(1/(1-theta)) ;
