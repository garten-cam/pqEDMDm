function dx = cstr_ode(~,x,d)
% cstr_ode returns the vector field for the simulation of a simple biochemical reactor

% x is the state
% d is the dilution tare

% Some constants
mu_max = 0.4; % Specific grouth rate
km = 0.12; % Monod saturation constant
k1 = 0.4545; % inverse of inhibition constant
Y = 0.4; % Process Yield r1/r2

% Assume that the species input concentration
% is zero, and that the substrate concentration
% at the input is
x2f = 4.0;

mu = mu_max * x(2) / (km + x(2) + k1 * x(2)^2);

dx = [
	x(1) * (mu - d);
	(x2f - x(2)) * d - (mu * x(1)) / Y
	];
end
