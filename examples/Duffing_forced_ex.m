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
figpath = "./pqEDMDm/examples/figures/";
% Example to work with forcing sgnals
rng(1) % For consistency
%%
% define the parameters for the simulation
num_ics = 9; % Number of initial conditions for the test
ics_width = 3; % ics range width
% Create the initial conditions for the orbits
ics = ics_width*rand(num_ics,2) - ics_width/2;
tfin = 20;
n_points = 601;
%%
% For the forced example I am only going to implement the duffing equation
% with multistability
% Two asymptotically stable points response
% parameters
tas.alpha = -1;
tas.beta = 1;
tas.delta = 0.5;
% The pqEDMD class accepts a structire array where the only necessary field
% in the state variables. It is not a tensor, because not all the
% trajectories are of the same lenght.

% preallocate the structure of tas orbits
tas_o = repmat(struct('y', zeros(n_points + 1, 2), ...
  'u', zeros(n_points + 1, 2),... forcing signals
  't', zeros(n_points + 1, 1)), num_ics,1);
% I am saving the 't' time array only for plotting purposes. The algorithm
% does not mind if that field is in there

% For the input, just a number, a frequency scaler
in = 10*rand(num_ics, 1) + 3;
% in = 2*rand(num_ics, 1) - 1;
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
  [tas_o(orb).t, tas_o(orb).y] = ode23(@(t,x)DuffEqODEu(t,x,tas, ...
    in(orb)),...
    0:tfin/n_points:tfin, ...
    ics(orb,:), ...
    odeSettings);
  % tas_o(orb).u = in(orb)*ones(size(tas_o(orb).t));%*cos(tas_o(orb).t);
  tas_o(orb).u = cos(in(orb)*tas_o(orb).t);
end

% Normalization 
% range = [-1,1];
% tas_n = dataset_normalization(tas_o,range);


%%
% Test the orthogonal pqEDMD
ts = [1 2]; % index of training trajectories
tr = 1:num_ics;
tr(ts) = [];
% create the decomposition object
tas_pq = pqEDMDm(p=[4], ... [3 4 5] 
               q=[1], ... [0.5 1 2]
               observable = @legendreObservable, ...
               dyn_dcp = @sidDecomposition); % '' to use the ordinary least squares
tas_ols = tas_pq.fit(tas_o(tr));
% The new iteration of the algorithm does not need a tr_ts thing. Just feed
% the ncessary training trajectories into the new fit function
% 
%
%%
% test the prdiction with the first sample
% load("x0.mat")
% test_one = tas_ols(1).predict(x0',size(tas_n(1).t,1)-11,{tas_n(1).u(11:end,:)})
% figure
% hold on
% plot(tas_n(1).sv(11:end,1),tas_n(1).sv(11:end,2),'b')
% plot(test_one.y(:,1),test_one.y(:,2),'k')
%%
% errors
% preallocate
err = zeros(numel(tas_ols),1);
for decp = 1 : numel(tas_ols)
  err(decp) = tas_ols(decp).error(tas_o(ts));
end
% where is the min?
[~, best] = min(err);
%%
% best = 4;
tas_p = tas_ols(best).pred_from_test(tas_o(ts));
%
% Plot
tas_f = figure(1);
clf
lay_tas = tiledlayout(3,3,"TileSpacing","tight");
for tr_i = 1 : numel(tr)
    nexttile(tr(tr_i))
    plot(tas_o(tr(tr_i)).t, tas_o(tr(tr_i)).y, 'b')
end
for ts_i = 1 : numel(ts)
    nexttile(ts(ts_i))
    hold on
    plot(tas_o(ts(ts_i)).t, tas_o(ts(ts_i)).y, 'r')
    plot(tas_o(ts(ts_i)).t, tas_p(ts_i).y, '-.k')
end
xlabel(lay_tas,'t','interpreter','latex')
ylabel(lay_tas,'$x_1$,$x_2$','interpreter','latex')
% saveas(tas_f,strcat(figpath, "forced_duff.png")),
% eigA = figure(2);
% clf
% hold on
% plot(cos(0:pi/50:2*pi),sin(0:pi/50:2*pi))
% scatter(real(eig(tas_ols(best).A)),imag(eig(tas_ols(best).A)),"filled")
% Even though the convergence is not perfect, this apprximation should be
% enough for control.


%%% Differential equation to solve
function Dx = DuffEqODEu(t,X,P,u)
%DuffEqODE
Dx1 = X(2);
Dx2 = -P.delta*X(2) - P.alpha*X(1) - P.beta*X(1)^3 + cos(u*t);
Dx = [Dx1;Dx2];
end

function data_n = dataset_normalization(data,range)
%  It is necessary to normalize. Get the normalization according to the
%  complete dataset
[~,y.ctr,y.scl] = normalize(cell2mat({data.y}'),"range",[-1,1]);
% the same for the input
[~,u.ctr,u.scl] = normalize(cell2mat({data.u}'),"range",range);
% apply the normalization to every sample. according t the center and the
% scale
data_n = arrayfun(@(x) struct('y', normalize(x.y,"center",y.ctr, ...
  "scale",y.scl), ...
  'u', normalize(x.u,"center",u.ctr, ...
  "scale",u.scl), ...
  't', x.t),data); % keep the time as is

end