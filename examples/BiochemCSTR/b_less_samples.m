% Script for continuing the simulation of the bioreactor.
% This time, we will have more "realistic" initial conditions
rng(4321)
% Define the parameters for the simulation
num_ics = 10;
% For this example, the initial conditions will be in the lower left corner. This exemplifies a low species and substrate concentration at the beginning of the experiment
ics = [(0.1-0.05)*rand(num_ics, 1)+0.05, (0.2-0.01)*rand(num_ics, 1)+0.01];
%
tfin = 40;
n_points = 1*tfin+1;
%
in = (0.16-0.1)*rand(num_ics, 1)+0.1;
% preallocate the array of structures
exp = arrayfun(@(z)struct('y', zeros(n_points, 2), ...
	'u', zeros(n_points, 1),... forcing signals
	't', zeros(n_points, 1)),1:num_ics)';
% Run the ODE
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[exp(orb).t, exp(orb).y] = ode23s(@(t,x)cstr_ode(t,x, ...
		in(orb)),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	exp(orb).u = in(orb)*ones(size(exp(orb).t));
end
% Normalixe the experiments
exp_n = normalize_data(exp,[-1,1]);
% Define the same wrapper. As the other script
pen_edmd = pqEDMDm( ...
	p = [2 4 5 6],...
	q = [0.5 1 1.5 2.0 2.5],...
	observable = @legendreObservable,...
	dyn_dcp = @svdDecomposition...
	);
% set the training and testing sets
tr = [1 2 6 7 8 9];
ts = [3 4 5 10];
% tr(ts) = [];
% Calculate all the decompositions
pen_dcps = pen_edmd.fit(exp_n(tr));
%
% Test the decompositions, calculate the error on the test set
err = arrayfun(@(dcp)dcp.error(exp_n(ts)),pen_dcps);
% Where is the min
[minerr, best] = min(err);
% extract the best
pen_dcp = pen_dcps(best);
% Predict from the test

pen_pred = pen_dcp.pred_from_test(exp_n(ts));

pq_fig = figure(4);
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
%%
% What happens if we add noise?
% Define a process noise, we do not know for sure the specific growth rate
% mu_max = .4 + .4/20 * rand(num_ics, 1);
% mu_max = normrnd(0.4,0.0,num_ics,1);
% preallocate the stochastic experiments
exp_stch = arrayfun(@(z)struct('y', zeros(n_points, 2), ...
	'u', zeros(n_points, 1),... forcing signals
	't', zeros(n_points, 1)),1:num_ics)';
% Run the simulation
for orb = 1 : num_ics
	[exp_stch(orb).t, exp_stch(orb).y_cl] = ode23s(@(t,x)cstr_ode(t,x, ...
		in(orb)),... there is a different mu_max per experiment
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	exp_stch(orb).y = exp_stch(orb).y_cl + normrnd(0,.01,size(exp_stch(orb).y));
	exp_stch(orb).u = in(orb)*ones(size(exp(orb).t));
end
% Normalize the experiments
exp_nstch = normalize_data(exp_stch,[-1,1]);
pen_edmd_stch = pqEDMDm( ...
	p = [2 4 5 6],...
	q = [0.5 1 1.5 2.0 2.5],...
	observable = @legendreObservable,...
	dyn_dcp = @svdDecomposition...
	);
% With the same svd decomposition, the approximation becomes
stch_dcps = pen_edmd_stch.fit(exp_nstch(tr));
% Calculate the error
st_err = arrayfun(@(dcp)dcp.error(exp_nstch(ts)),stch_dcps);
% whos the best
[st_min, st_bt] = min(st_err);
% extract the best
st_dcp = stch_dcps(st_bt);
% Make the prediction from the test set
stch_pred = st_dcp.pred_from_test(exp_nstch(ts));


% Plot the result
st_fig = figure(5);
% instead of a phase plane, use a grid
clf
hold on
st_trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_nstch(tr));
% Plot the testing set
st_tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_nstch(ts));
% Plot the testing approximation
tsa = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),stch_pred);
xlabel('$x_1$',Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')
legend([trp(1), tsp(1), tsa(1)],{'training','testing','pqEDMD'},"Location","best")
title("pqEDMD appx")
%%
% It is nice, but, can we do better?
% Create a pqEDMD with sidDecomposition
% sid_edmd = pqEDMDm( ...
% 	p = [2 3 4 5],...
% 	q = [0.5 1 1.5 2.0 2.5],...
% 	observable = @legendreObservable,...
% 	dyn_dcp = @sidDecomposition...
% 	);
% % With the new sid decomposition, the approximation becomes
% sid_dcps = sid_edmd.fit(exp_stch(tr));
% % Calculate the error
% sid_err = arrayfun(@(dcp)dcp.error(exp_stch(ts)),sid_dcps);
% % whos the best
% [sid_min, sid_bt] = min(sid_err);
% % extract the best
% sid_dcp = sid_dcps(sid_bt);
% % Make the prediction from the test set
% sid_pred = sid_dcp.pred_from_test(exp_stch(ts));
% % Plot
% sid_fig = figure(6);
% % instead of a phase plane, use a grid
% clf
% hold on
% sid_trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_stch(tr));
% % Plot the testing set
% sid_tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_stch(ts));
% % Plot the testing approximation
% sid_ts = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),sid_pred);
% xlabel('$x_1$',Interpreter='latex')
% ylabel('$x_2$', Interpreter='latex')
% legend([trp(1), tsp(1), sid_ts(1)],{'training','testing','pqEDMD'},"Location","best")
% title("sidEDMD appx")
