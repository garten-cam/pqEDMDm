function dx = cstr_ode_stch(~, x, d, mu_max)
% cstr_ode_stch returns the vector field of the stochastic
% bioreactor. For this case, the dynamics are unknown
%
% x is the state
% d is the dilution rate
%
% define the constants
km = 0.12; % Monod saturation constant
k1 = 0.4545; % inverse of inhibition constant
Y = 0.4; % Process Yield r1/r2
% Again, this assumes that there is no concentration of species in the input
x2f = 4.0;
% Calculate mu this is the only difference with the deterministic case,
% here, mu_max is random variable.
mu = mu_max * x(2) / (km + x(2) + k1 * x(2)^2);

dx = [x(1) * (mu - d); (x2f - x(2)) * d - (mu * x(1)) / Y];
end

