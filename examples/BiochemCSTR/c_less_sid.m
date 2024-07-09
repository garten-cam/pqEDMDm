% After cutting on the amount of samples, and restricting the initial conditions
% to the lower left corner, lets try the sid decomposition on the same set.
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

% Define the same wrapper. As the other script
sid_edmd = pqEDMDm( ...
	p = [2 3 4],...
	q = [0.5 0.8 1 1.5],...
	observable = @legendreObservable,...
	dyn_dcp = @sidDecomposition...
	);
% set the training and testing sets
tr = [1 2 6 7 8 9];
% tr = [1 2];
ts = [3 4 5 10];
% fit it to the training set
sid_dcps = sid_edmd.fit(exp_n(tr));
% errors for test set
err_sid = arrayfun(@(dcp)dcp.error(exp_n(ts)), sid_dcps);
% who is the best?
[~, best_sid] = min(err_sid);
% Extract the best
sid_dcp = sid_dcps(best_sid);
% Make the prediction on the test set
sid_pred = sid_dcp.pred_from_test(exp_n(ts));
% plot the result
sid_fig = figure(6);
clf
hold on
% training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_n(tr));
% testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_n(ts));
% and the approxx
sip = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),sid_pred);

