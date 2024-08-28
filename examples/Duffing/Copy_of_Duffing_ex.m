figpath = "./pqEDMDm/examples/figures/";
% for numerical consistency
rng(1)
% define the parameters for the simulation
num_ics = 10; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbitst
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 30;
n_points = 10*tfin + 1;
% Two asymptotically stable points response
% parameters
tas.alpha = -1;
tas.beta = 1;
tas.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.
% %
% preallocate the structure of tas orbits
two_as = repmat(struct('y', zeros(n_points, 2), ...
	't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[two_as(orb).t, two_as(orb).y] = ode23(@(t,x)DuffEqODE(t,x,tas),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
end

two_nrm = normalize_data(two_as,[-1,1]);
% two_nrm = two_as;
  


%%
% Test the orthogonal pqEDMD
ts = [1 2]; % index of training trajectories
tr = [3 4 5 6 7 8 9 10];
% create the decomposition object
tas_pq = pqEDMDm(p=[2 3 4 5], ...
	q=[0.5 1 1.5 2], ...
	observable = @legendreObservable, ...
	dyn_dcp = @svdDecomposition); % '' to use the ordinary least squares
dcps = tas_pq.fit(two_nrm(tr));
%%
% Get the best performing decomposition
err = zeros(numel(dcps),1);
for decp = 1 : numel(dcps)
	err(decp) = dcps(decp).abs_error(two_nrm(ts));
end

[~, best] = min(err);
dcp = dcps(best);

% Compare the prediction of the best performing prediction with
% the test set 


appx = dcp.pred_from_test(two_nrm(ts));

%
tas_f = figure(1);
clf
hold on
% plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b',LineWidth=2),two_nrm(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r',LineWidth=2),two_nrm(ts));
% Plot the sid approximation
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=2),appx);
% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))