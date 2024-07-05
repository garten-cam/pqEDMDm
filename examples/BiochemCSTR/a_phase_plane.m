% Script for the simulation of the simple bioreactor
% simulate the phase plane and the approximation by the pqEDMD
rng(12345)
% define the parameters for the simulation
num_ics = 60; % Number of initial conditions for the test
% Create the initial conditions for the samples
% Even though the initial conditions should always start at the lower left
% corner, this will plot the phase plane with constant dilution rate.
% (b-a)*rand(n,m)+a
% For $x_0\in (0,\;1.58]\times(0,\;4]$
ics = [
	(1.58-0.01)*rand(num_ics, 1)+0.01, (4-0.01)*rand(num_ics, 1)+0.01
	];
%
tfin = 40;
n_points = 1*tfin+1;
%
% The first figure is the phase plane with dilution rate D=O.2
in = 0.2*ones(num_ics,1);
% Preallocate the experiments if the samples are not in a column
% vector, the algorithm breaks. Dammit
exp = arrayfun(@(z)struct('y', zeros(n_points, 2), ...
	'u', zeros(n_points, 1),... forcing signals
	't', zeros(n_points, 1)),1:num_ics)';


odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[exp(orb).t, exp(orb).y] = ode23s(@(t,x)cstr_ode(t,x, ...
		in(orb)),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	exp(orb).u = in(orb)*ones(size(exp(orb).t));
end

%
trfig = figure(1);
clf
hold on
arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b'),exp)
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
title("Phase plane")

% Create the decomposition object
pen_edmd = pqEDMDm( ...
	p = [4 5 6],...
	q = [0.5 1 2 2.5],...
	observable = @legendreObservable,...
	dyn_dcp = @svdDecomposition...
	);
% Normalize the data
exp_n = normalize_data(exp,[-1,1]);
% set the training and testing sets
tr = 1:2:num_ics;
ts = 2:2:num_ics;
% tr(ts) = [];
% Calculate all the decompositions
pen_dcps = pen_edmd.fit(exp_n(tr));
%
% Test the decompositions, calculate the error on the test set
err = arrayfun(@(dcp)dcp.error(exp_n(ts)),pen_dcps);
% where is the min?
[minerr, best] = min(err);
% extract the best
pen_dcp = pen_dcps(best);
% Predict from the test

pen_pred = pen_dcp.pred_from_test(exp_n(ts));

pq_fig = figure(2);
clf
hold on
% Plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_n(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_n(ts));
% Plot the testing approximation
tsa = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),pen_pred);
xlabel('$x_1$',Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')
legend([trp(1), tsp(1), tsa(1)],{'training','testing','pqEDMD'},"Location","best")
title("pqEDMD appx")

% Put the two figures in a tiled layout for the paper
a_pp = figure(3);
clf
tiledlayout("horizontal","TileSpacing","tight")
nexttile
hold on
arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b'),exp)
xlabel('$x_1$: biomass', 'Interpreter', 'latex')
ylabel('$x_2$: substrate', 'Interpreter', 'latex')
title("(a) Phase plane")
nexttile
hold on
% Plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_n(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_n(ts));
% Plot the testing approximation
tsa = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),pen_pred);
xlabel('$x_1$: biomass',Interpreter='latex')
ylabel('$x_2$ substrate', Interpreter='latex')
legend([trp(1), tsp(1), tsa(1)],{'training','testing','pqEDMD'})
title("(b) pqEDMD appx")
set(gcf,'PaperPosition', [0, 0, 8, 4])
% saveas(a_pp,"~/Documents/sidDecompositionJ/figures/a_phase_plane.fig")
% saveas(a_pp,"~/Documents/sidDecompositionJ/figures/a_phase_plane.eps", "epsc")
