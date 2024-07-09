% figpath = "./pqEDMDm/examples/figures/";
% Example to work with forcing sgnals
rng(1) % For consistency
%%
% define the parameters for the simulation
num_ics = 12; % Number of initial conditions for the test
ics_width = 2; % ics range width
% Create the initial conditions for the orbits
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 20;
n_points = 20*tfin+1;
%%
% For the forced example I am only going to implement the duffing equation
% with multistability
% Two asymptotically stable points response
% parameters
duffor.alpha = -1;
duffor.beta = 1;
duffor.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of the experiments
duff_exp = arrayfun(@(z)struct('y', zeros(n_points, 2), ...
	'u', zeros(n_points, 2),... forcing signals
	't', zeros(n_points, 1)), 1:num_ics)';

% For the input, just a number, a frequency scaler
omega_in = 10*rand(num_ics, 1) + 3;
in = arrayfun(@(om_in)struct('omega',om_in,'gamma',1),omega_in);
% in = 2*rand(num_ics, 1) - 1;
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[duff_exp(orb).t, duff_exp(orb).y] = ode23(@(t,x)DuffEqODEu(t,x,duffor, ...
		in(orb)),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	duff_exp(orb).u = in(orb).gamma*cos(in(orb).omega*duff_exp(orb).t);
end

% Normalization
range = [-1,1];
exp_n = normalize_data(duff_exp,range);


%%
ts = [2 7]; % index of training trajectories
tr = 1:num_ics;
tr(ts) = [];
% create the decomposition object
tas_pq = pqEDMDm(p=[2 3 4], ... [3 4 5]
	q=[0.5 1 1.3], ... [0.5 1 2]
	observable = @legendreObservable, ...
	dyn_dcp = @sidDecomposition); % '' to use the ordinary least squares
% Fit the decompositions
duff_dcps = tas_pq.fit(exp_n(tr)); % All decompositions
%
% errors
err = arrayfun(@(dcp)dcp.error(exp_n(ts)), duff_dcps);
% where is the min?
[~, best] = min(err);
%
% Extract the best
duff_dcp = duff_dcps(best);
duff_pred = duff_dcp.pred_from_test(duff_exp(ts));
%
% Plot
tas_f = figure(1);
clf
lay_tas = tiledlayout(4,3,"TileSpacing","tight");
for tr_i = 1 : numel(tr)
	nexttile(tr(tr_i))
	plot(exp_n(tr(tr_i)).t, exp_n(tr(tr_i)).y, 'b')
end
for ts_i = 1 : numel(ts)
	nexttile(ts(ts_i))
	hold on
	plot(exp_n(ts(ts_i)).t, exp_n(ts(ts_i)).y, 'r')
	plot(exp_n(ts(ts_i)).t, duff_pred(ts_i).y, '-.k')
end
xlabel(lay_tas,'t','interpreter','latex')
ylabel(lay_tas,'$x_1$,$x_2$','interpreter','latex')
% saveas(tas_f,strcat(figpath, "forced_duff.png")),
sid_fig = figure(2);
clf
hold on
% training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b','LineWidth',2),exp_n(tr));
% testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r','LineWidth',2),exp_n(ts));
% and the approxx
sip = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k','LineWidth',2),duff_pred);
