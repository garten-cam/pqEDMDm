% Script to generate the samples for testing the different algorithms using
% the unforced and stochastic Duffing equation.
% The stochastic part is in the measurement

% \begin{align}
% 	\dot{x}_1(t)&=x_2(t)\\
% 	\dot{x}_2(t)&=-\delta{}x_2(t)-x_1(t)(\beta+\alpha{}x_1^{2}(t)). \\
% y(t) &=x(t) + w(t)
% \end{align}

% In this script we will test the case that has two assymptotically stable
% points. Then,
% \delta = 0.5
% \alpha = -1
% \beta = 1
%-------------------------------------------------------------------------%
% figpath = "./pqEDMDm/examples/figures/";
% for numerical consistency
rng(1)
% define the parameters for the simulation
num_ics = 6; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbitst
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 40;
n_points = 10*tfin + 1;
% Two asymptotically stable points response
% parameters
p.alpha = -1;
p.beta = 1;
p.delta = 0.5;
% The noise
sigma = 0.1;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.
% %
% preallocate the structure of the orbits or experiments
duff_exp = repmat(struct('y', zeros(n_points, 2), ...
	't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there
% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[duff_exp(orb).t, y] = ode23(@(t,x)DuffEqODE(t,x,p),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	duff_exp(orb).y_det = y;
	duff_exp(orb).y = y + normrnd(0, sigma, size(y));
end

% From here, we can call the  different decompositions
