% Script to generate the samples for testing the different algorithms using
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
% figpath = "./pqEDMDm/examples/figures/";
% for numerical consistency
rng(1)
% define the parameters for the simulation
num_ics = 9; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbitst
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 30;
n_points = 8*tfin + 1;
% Two asymptotically stable points response
% parameters
p.alpha = -1;
p.beta = 1;
p.delta = 0.5;
% Input parameters. In DuffEqODEu, the input is u(t)=gamma*cos(omega*t)
u = arrayfun(@(gammai)struct('gamma',gammai,'omega',0),4*rand(num_ics,1)-2);
% preallocate the output
duff_exp = arrayfun(@(exp)struct('t',zeros(n_points,1),...
  'y',zeros(n_points,2),'u',exp.gamma*ones(n_points,1)),u);
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[duff_exp(orb).t, y] = ode23(@(t,x)DuffEqODEu(t,x,p, ...
		u(orb)),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...  
		odeSettings);
  duff_exp(orb).y_det = y;
  duff_exp(orb).y = y + normrnd(0, 0.02, size(y));
end

% From here, call different decompositions to test the performance
