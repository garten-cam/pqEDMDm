%%
% The last tutrorial with the orbits class was very confusing. This new one
% should be more clear.
%%
rng(1)
%%
% define the parameters for the simulation
num_ics = 6; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbits
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 20;
n_points = 10*tfin;
%%
% Undamped response
% parameters
und.alpha = 1;
und.beta = 1;
und.delta = 0;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of undamped orbits
und_o = repmat(struct('sv', zeros(n_points, 2), ...
                            't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [und_o(orb).t, und_o(orb).sv] = ode23(@(t,x)DuffEqODE(t,x,und),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);       
end
%%
% Test the orthogonal pqEDMD
tr = 3; % index of training trajectories
ts = [1 2 4 5 6];
% create the decomposition object
und_pq = pqEDMD(p=[2 3 4 5], ...
               q=[0.5 1 2], ...
               polynomial='Legendre', ...
               method=''); % '' to use the ordinary least squares
und_ols = und_pq.fit(und_o(tr));
% The new iteration of the algorithm does not need a tr_ts thing. Just feed
% the ncessary training trajectories into the new fit function
%
%%
% errors
% preallocate
err = zeros(numel(und_ols),1);
for decp = 1 : numel(und_ols)
err(decp) = und_ols(decp).error(und_o(ts));
end
% where is the min?
[~, best] = min(err);
% Get the predictions from the EDMD models
%%
% !!preallocate all the things!!!
% array of points per test trajectory (they could be different). In this
% case they are the same.
num_p = arrayfun(@(x) size(x.sv,1), und_o(ts));
% loop and get the predictions

und_p = und_ols(best).predict(ics(ts,:), num_p);
%%
% Plot the result
und_f = figure(1);
clf
hold on
for tr_i = 1 : numel(tr)
    plot(und_o(tr(tr_i)).sv(:,1), und_o(tr(tr_i)).sv(:,2), 'r')
end
for ts_i = 1 : numel(ts)
    plot(und_o(ts(ts_i)).sv(:,1), und_o(ts(ts_i)).sv(:,2), 'b')
    plot(und_p(ts_i).sv(:,1),und_p(ts_i).sv(:,2), '--k')
end

%% -------------------------------------------------- %%
% Galobaly asymptotically stable response
% parameters
gas.alpha = 1;
gas.beta = 1;
gas.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of gas orbits
gas_o = repmat(struct('sv', zeros(n_points, 2), ...
                            't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [gas_o(orb).t, gas_o(orb).sv] = ode23(@(t,x)DuffEqODE(t,x,gas),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);       
end
%
% Test the orthogonal pqEDMD
tr = 5; % index of training trajectories
ts = [1 2 3 4 6];
% create the decomposition object
gas_pq = pqEDMD(p=[2 3 4 5], ...
               q=[0.5 1 2], ...
               polynomial='Legendre', ...
               method=''); % '' to use the ordinary least squares
gas_ols = gas_pq.fit(gas_o(tr));
% The new iteration of the algorithm does not need a tr_ts thing. Just feed
% the ncessary training trajectories into the new fit function
%
%
% errors
% preallocate
err = zeros(numel(gas_ols),1);
for decp = 1 : numel(gas_ols)
err(decp) = gas_ols(decp).error(gas_o(ts));
end
% where is the min?
[~, best] = min(err);
% Get the predictions from the EDMD models
%
% !!preallocate all the things!!!
% array of points per test trajectory (they could be different). In this
% case they are the same.
num_p = arrayfun(@(x) size(x.sv,1), gas_o(ts));
% loop and get the predictions
gas_p = gas_ols(best).predict(ics(ts,:), num_p);
%
% Plot the result
gas_f = figure(2);
clf
hold on
for tr_i = 1 : numel(tr)
    plot(gas_o(tr(tr_i)).sv(:,1), gas_o(tr(tr_i)).sv(:,2), 'r')
end
for ts_i = 1 : numel(ts)
    plot(gas_o(ts(ts_i)).sv(:,1), gas_o(ts(ts_i)).sv(:,2), 'b')
    plot(gas_p(ts_i).sv(:,1),gas_p(ts_i).sv(:,2), '--k')
end

%% -------------------------------------------------- %%
% Two asymptotically stable points response
% parameters
tas.alpha = -1;
tas.beta = 1;
tas.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of tas orbits
tas_o = repmat(struct('sv', zeros(n_points, 2), ...
                            't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Loop for all initial conditions
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [tas_o(orb).t, tas_o(orb).sv] = ode23(@(t,x)DuffEqODE(t,x,tas),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);       
end
%
% Test the orthogonal pqEDMD
tr = 3; % index of training trajectories
ts = [1 2 4 5 6];
% create the decomposition object
tas_pq = pqEDMD(p=[2 3 4 5], ...
               q=[0.5 1 2], ...
               polynomial='Legendre', ...
               method=''); % '' to use the ordinary least squares
tas_ols = tas_pq.fit(tas_o(tr));
% The new iteration of the algorithm does not need a tr_ts thing. Just feed
% the ncessary training trajectories into the new fit function
%
%
% errors
% preallocate
err = zeros(numel(tas_ols),1);
for decp = 1 : numel(tas_ols)
err(decp) = tas_ols(decp).error(tas_o(ts));
end
% where is the min?
[~, best] = min(err);
% Get the predictions from the EDMD models
%
% !!preallocate all the things!!!
% array of points per test trajectory (they could be different). In this
% case they are the same.
num_p = arrayfun(@(x) size(x.sv,1), tas_o(ts));
% loop and get the predictions
tas_p = tas_ols(best).predict(ics(ts,:), num_p);
%
% Plot the result
tas_f = figure(3);
clf
hold on
for tr_i = 1 : numel(tr)
    plot(tas_o(tr(tr_i)).sv(:,1), tas_o(tr(tr_i)).sv(:,2), 'r')
end
for ts_i = 1 : numel(ts)
    plot(tas_o(ts(ts_i)).sv(:,1), tas_o(ts(ts_i)).sv(:,2), 'b')
    plot(tas_p(ts_i).sv(:,1),tas_p(ts_i).sv(:,2), '--k')
end
% !! Taking the outermost trajectory of the set, it is able to generalize
% the other side!! 



























% Differential equation to solve
function Dx = DuffEqODE(~,X,P)
%DuffEqODE
Dx1 = X(2);
Dx2 = -P.delta*X(2) - P.alpha*X(1) - P.beta*X(1)^3;
Dx = [Dx1;Dx2];
end