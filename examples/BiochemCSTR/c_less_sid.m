% After cutting on the amount of samples, and restricting the initial conditions
% to the lower left corner, lets try the sid decomposition on the same set.
clear variables
rng(4321)
% Define the parameters for the simulation
num_ics = 10;
% For this example, the initial conditions will be in the lower left corner. This exemplifies a low species and substrate concentration at the beginning of the experiment
ics = [(0.1-0.05)*rand(num_ics, 1)+0.05, (0.2-0.01)*rand(num_ics, 1)+0.01];
%
tfin = 40;
n_points = 1*tfin+1;
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

% Normalize the data
exp_n = normalize_data(exp,[-1,1]);
% exp_n = exp;
% Define the same wrapper. As the other script
sid_edmd = pqEDMDm( ...
	p = [2 3 4 5],...
	q = [0.2 0.5 1 1.5 2.0 2.5],...
	observable = @hermiteObservable,...
	dyn_dcp = @(obs,sys)sidDecomposition(4,1,obs,sys)...
	);
% set the training and testing sets
tr = [1 2 6 7 8 9];
ts = [3 4 5 10];
% fit it to the training set
sid_dcps = sid_edmd.fit(exp_n(tr));
% errors for test set
err_sid = arrayfun(@(dcp)dcp.error(exp_n(ts)), sid_dcps);
% who is the best?
[~, best_sid] = min(err_sid);
% Extract the best
%
sid_dcp = sid_dcps(best_sid);
% Make the prediction on the test set
sid_pred = sid_dcp.pred_from_test(exp_n(ts));
% plot the result
sid_fig = figure(7);
clf
hold on
% training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_n(tr));
% testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_n(ts));
% and the approxx
sip = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),sid_pred);
%%
% Now with  noise
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
	exp_stch(orb).y = exp_stch(orb).y_cl + normrnd(0,.005,size(exp_stch(orb).y));
	exp_stch(orb).u = in(orb)*ones(size(exp(orb).t));
end
% Normalize the experiments
%%
exp_nstch = normalize_data(exp_stch,[-1,1]);
pen_edmd_stch = pqEDMDm( ...
  p = [3 4],...
	q = [0.5 1 1.5],...
	observable = @legendreObservable,...
	dyn_dcp = @(s,o)sidDecomposition(1,2,s,o)...
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
st_fig = figure(8);
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