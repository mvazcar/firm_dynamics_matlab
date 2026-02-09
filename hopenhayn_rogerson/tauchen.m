% This function discretizes the continuous AR(1) process
% by using the method proposed by Tauchen (1986). 
%The AR(1) process takes the following form:
% y(t) = (1-rho)*mu + rho*y(t-1) + eps(t), where eps ~ N(0,sig^2).

function [s, Pi] = tauchen(mu,rho,sig,N)

m       = 5 ;
s       = zeros(N,1);
Pi      = zeros(N,N);

s(1)    = mu - m*sqrt(sig^2/(1-rho^2));
s(N)    = mu + m*sqrt(sig^2/(1-rho^2));

step    = (s(N)-s(1))/(N-1);

for i=2:(N-1)
   s(i) = s(i-1) + step;
end

for j = 1:N
    for k = 1:N
        if k == 1
            if ~isreal((s(1) - (1-rho)*mu - rho*s(j) + step/2) / sig)
            keyboard
            end
            Pi(j,k) = cdf_normal((s(1) - (1-rho)*mu - rho*s(j) + step/2) / sig) ;
        elseif k == N
            Pi(j,k) = 1 - cdf_normal((s(N) - (1-rho)*mu - rho*s(j) - step/2) / sig);
        else
            Pi(j,k) = cdf_normal((s(k) - (1-rho)*mu - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - (1-rho)*mu - rho*s(j) - step/2) / sig);
        end
    end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));
