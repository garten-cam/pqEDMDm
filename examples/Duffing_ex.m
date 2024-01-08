figpath = "./pqEDMDm/examples/figures/";
% for numerical consistency
rng(1)
% define the parameters for the simulation
num_ics = 6; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbits
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 20;
n_points = 10*tfin;
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
tas_o = repmat(struct('y', zeros(n_points, 2), ...
                      't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [tas_o(orb).t, tas_o(orb).y] = ode23(@(t,x)DuffEqODE(t,x,tas),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);       
end



%%
% Test the orthogonal pqEDMD
tr = [3 5]; % index of training trajectories
ts = [1 2 4 6];
% create the decomposition object
tas_pq = pqEDMDm(p=[5 6 7 ], ...
               q=[0.5 1 2], ...
               observable = @legendreObservable, ...
               dyn_dcp = @pqDecomposition); % '' to use the ordinary least squares
tas_ols = tas_pq.fit(tas_o(tr));
%%
% Get the best performing decomposition
err = zeros(numel(tas_ols),1);
for decp = 1 : numel(tas_ols)
err(decp) = tas_ols(decp).error(tas_o(ts));
end

[~, best] = min(err);
best_tas_ols = tas_ols(best);

% Compare the prediction of the best performing prediction with 
% the test set


tas_approx = best_tas_ols.pred_from_test(tas_o(ts));


tas_f = figure(1);
clf
hold on
for tr_i = 1 : numel(tr)
    plot(tas_o(tr(tr_i)).y(:,1), tas_o(tr(tr_i)).y(:,2), 'b')
end
for ts_i = 1 : numel(ts)
    plot(tas_o(ts(ts_i)).y(:,1), tas_o(ts(ts_i)).y(:,2), 'r')
end
axis([-2.1 2.1 -2.1 2.1])
saveas(tas_f,strcat(figpath, "tr_ts.png"))
for ts_i = 1 : numel(ts)
    plot(tas_approx(ts_i).y(:,1),tas_approx(ts_i).y(:,2), '.k')
end
saveas(tas_f,strcat(figpath, "approx.png"))

% Differential equation to solve
function Dx = DuffEqODE(~,X,P)
%DuffEqODE
Dx1 = X(2);
Dx2 = -P.delta*X(2) - P.alpha*X(1) - P.beta*X(1)^3;
Dx = [Dx1;Dx2];
end