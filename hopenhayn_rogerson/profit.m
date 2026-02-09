% ##########################################################################
% Optimal Firm Profit - Static
% ##########################################################################

function f = profit(params, s, p, w)

theta = params.theta ;
cf = params.cf ;

prof = (p*exp(s)*theta^theta).^(1/(1-theta))*(1 - theta) - w*cf ;
f = prof ;