% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2023 Camilo Garcia-Tenorio
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this
% software and associated documentation files (the "Software"), to deal in the Software
% without restriction, including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or
% substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.

% Example to work with forcing sgnals and kalman filter
% run the forced example 
% define the parameters for the simulation
% rng(1)
num_ics = 6; % Number of initial conditions for the test
ics_width = 4; % ics range width
% Create the initial conditions for the orbits
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 20;
% the sampling frequency is 5 times faster than the forced_example
n_points = 50*tfin; 
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
    'u', zeros(n_points, 2),... forcing signals
    't', zeros(n_points, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% Create the input signal. for simplicity, it will be a step response,
% i.e., a constant input.
in = 100*rand(num_ics, 1) - 50;
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
    [tas_o(orb).t, tas_o(orb).sv] = ode23(@(t,x)DuffEqODEu(t,x,tas, ...
        in(orb)),...
        0:tfin/n_points:tfin, ...
        ics(orb,:), ...
        odeSettings);    
    tas_o(orb).u = ones(n_points + 1, 1)*in(orb);
end

%%
% Test the orthogonal pqEDMD
tr = [3 6]; % index of training trajectories
ts = [1 2 4 5];
% create the decomposition object
tas_pq = pqEDMD(p=[7 8], ...
               q=[0.5 1 2], ...
               polynomial='Legendre', ...
               method='rrr'); % '' to use the ordinary least squares
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
%%
tas_p = tas_ols(best).pred_from_test(tas_o(ts));
%%
% Plot
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
% Even though the convergence is not perfect, this apprximation should be
% enough for control.


%%% Differential equation to solve
function Dx = DuffEqODEu(~,X,P,u)
%DuffEqODE
Dx1 = X(2);
Dx2 = -P.delta*X(2) - P.alpha*X(1) - P.beta*X(1)^3 + u;
Dx = [Dx1;Dx2];
end