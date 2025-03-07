function duff_exp = duff_det_forced_step(num_ics, tfin, point_multiplier)
% function to generate the samples for testing the different algorithms using
% the forced and deterministic Duffing equation.
% The stochastic part is in the measurement

% \begin{align}
% 	\dot{x}_1(t)&=x_2(t)\\
% 	\dot{x}_2(t)&=-\delta{}x_2(t)-x_1(t)(\beta+\alpha{}x_1^{2}(t)) + u(t). \\
% y(t) &=x(t)
% \end{align}

% In this script we will test the case that has two assymptotically stable
% points. Then,
% \delta = 0.5
% \alpha = -1
% \beta = 1
%-------------------------------------------------------------------------%
ics_width = 4; % ics range width
% Create the initial conditions for the orbitst
ics = ics_width*rand(num_ics,2) - ics_width/2;
% tfin = 20;
n_points = point_multiplier * tfin + 1;
% Two asymptotically stable points response
% parameters
p.alpha = -1;
p.beta = 1;
p.delta = 0.5;
% Input parameters. In DuffEqODEu, the input is u(t)=gamma*cos(omega*t)
u = arrayfun(@(gammai)struct('gamma',gammai,'omega',0),4*rand(num_ics,1)-2);
% preallocate the output
duff_exp = arrayfun(@(u_i, np_i)struct( ...
	't',zeros(np_i,1),...
  'y',zeros(np_i, 2), ...
	'u',u_i.gamma*ones(np_i,1)),u,n_points);
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[duff_exp(orb).t, duff_exp(orb).y] = ode23(@(t,x)DuffEqODEu(t,x,p, ...
		u(orb)),...
		linspace(0,tfin(orb),n_points(orb)), ...
		ics(orb,:), ...  
		odeSettings);
end

% From here, call different decompositions to test the performance
