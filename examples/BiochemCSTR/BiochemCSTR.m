% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2024 Camilo Garcia-Tenorio
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

% Script for the simulation of the simple biological reactor. 
rng(54321)
% define the parameters for the simulation
num_ics = 10; % Number of initial conditions for the test
% Create the initial conditions for the samples
% The first variable is the species. It should always start in a low value
% The second variable is the sybstrate. It shloud occupy most (alamost all)
% of the volume.
% (b-a)*rand(n,m)+a
ics = [
(0.5-0.01)*rand(num_ics, 1)+0.01, (1-0.01)*rand(num_ics, 1)+0.01
];
%
tfin = 50;
n_points = 10*tfin;
%
% Create the input signal. for simplicity, it will be a step response,
% i.e., a constant input in the range [0.1 - 1]
in = (0.28-0.08)*rand(num_ics, 1)+0.08;
% in = rand(num_ics, 1);

% in = 0.08*ones(num_ics, 1);

% Preallocate the experiments
exp = repmat(struct('y', zeros(n_points + 1, 2), ...
  'u', zeros(n_points + 1, 1),... forcing signals
  't', zeros(n_points + 1, 1)), num_ics,1);


odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
  [exp(orb).t, exp(orb).y] = ode23s(@(t,x)biochem(t,x, ...
    in(orb)),...
    0:tfin/n_points:tfin, ...
    ics(orb,:), ...
    odeSettings);
  exp(orb).u = in(orb)*ones(size(exp(orb).t));%*cos(cycles(orb).t);
  % cycles(orb).u = cos(in(orb)*cycles(orb).t);
end

% some plots to see how we are progressing.
% select now a training and a testing set
tr = [1 2 3 4 6 7 9 10];
ts = [5 8]; % The second argumet is the num of test cycles

% trfig = figure(1);
% clf
% tile_lay = tiledlayout(2,2,'TileSpacing','tight');
% for tri = 1 : numel(tr)
%     nexttile
%     plot(cycles(tr(tri)).t, cycles(tr(tri)).y)
% end

st_sp = figure(1);
clf
hold on
for cyc = tr
    trp = plot(exp(cyc).y(:,1), exp(cyc).y(:,2),'b',LineWidth=2);
end
for cyc = ts
    tsp = plot(exp(cyc).y(:,1), exp(cyc).y(:,2),'r',LineWidth=2);
end

xlabel('$x_1$',Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')
% For some reason, the first seed that I tested, gave 
% a set of orbits that match my criterion. Five and five.

% Create the decomposition object
pen_pq = pqEDMDm( ...
    p = [2 3 4 5],...
    q = [0.5 1 2 2.5],...
    observable = @legendreObservable,...
    dyn_dcp = @pqDecomposition...
    );

% Calculate the decomposition
pen_svd = pen_pq.fit(exp(tr));
%
% Test the decompositions
err = zeros(numel(pen_svd),1);
for decp = 1 : numel(pen_svd)
  err(decp) = pen_svd(decp).error(exp(ts));
end
% where is the min?
[minerr, best] = min(err)

% Predict from the test 
pen_p = pen_svd(best).pred_from_test(exp(ts));

% Plot the result
figure(1)
for cyc = 1 : numel(ts)
    tsa = plot(pen_p(cyc).y(:,1), pen_p(cyc).y(:,2),'-.k',LineWidth=2);
end
legend([trp, tsp, tsa],{'training','testing','EDMD'})
cleanfigure;
matlab2tikz('./Penicillin/figures/EDMDtrts.tex',...
    'width','5cm',...
    'height','5cm');

% Everything seems so nice...
% But wait, in "reality", systems are not
% deterministic. Lets try with some noice
%%
% Ater the simulation of the EDMD, test with
% noise.
mu_max = .4 + .4/10 * rand(num_ics, 1);
% 
exp_stch = repmat(struct('y', zeros(n_points + 1, 2), ...
  'u', zeros(n_points + 1, 1),... forcing signals
  't', zeros(n_points + 1, 1)), num_ics,1);


odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
  [exp_stch(orb).t, exp_stch(orb).y] = ode23s(@(t,x)biochem_stch(t,x, ...
      in(orb),...
      0.4),...
    0:tfin/n_points:tfin, ...
    ics(orb,:), ...
    odeSettings);
  exp_stch(orb).y = exp_stch(orb).y + exp_stch(orb).y/20 .* rand(size(exp_stch(orb).y));
  exp(orb).u = in(orb)*ones(size(exp(orb).t));%*cos(cycles(orb).t);
  % cycles(orb).u = cos(in(orb)*cycles(orb).t);
end

st_sp = figure(2);
clf
hold on
for cyc = tr
    trp = plot(exp_stch(cyc).y(:,1), exp_stch(cyc).y(:,2),'b',LineWidth=2);
end
for cyc = ts
    tsp = plot(exp_stch(cyc).y(:,1), exp_stch(cyc).y(:,2),'r',LineWidth=2);
end
xlabel('$x_1$',Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')

pen_pq_stch = pqEDMDm( ...
    p = [2 3 4 5],...
    q = [0.5 1 2 2.5],...
    observable = @legendreObservable,...
    dyn_dcp = @maxLikeDecomposition...
    );

% Calculate the decomposition
pen_ml = pen_pq_stch.fit(exp(tr));
%
% Test the decompositions
err_ml = zeros(numel(pen_ml),1);
for decp = 1 : numel(pen_ml)
  err_ml(decp) = pen_ml(decp).error(exp_stch(ts));
end
% where is the min?
[minerr_ml, best_ml] = min(err_ml)

% Predict from the test 
pen_p_stch = pen_ml(best).pred_from_test(exp_stch(ts));

% Plot the result
figure(2)
for cyc = 1 : numel(ts)
    tsa = plot(pen_p_stch(cyc).y(:,1), pen_p_stch(cyc).y(:,2),'-.k',LineWidth=2);
end
legend([trp, tsp, tsa],{'training','testing','EDMD'})
%%
function dx = biochem_stch(~, x, d, mu_max)
% x is the state
% d is the dilution tare

% Some constants
% mu_max = 0.4; % Specific grouth rate
km = 0.12; % Monod saturation constant
k1 = 0.4545; % inverse of inhibition constant
Y = 0.4; % Process Yield r1/r2

% Assume that the species input concentration 
% is zero, and that the substrate concentration
% at the input is
x2f = 4.0;

mu = mu_max * x(2) / (km + x(2) + k1 * x(2)^2);

dx = [
x(1) * (mu - d);
(x2f - x(2)) * d - (mu * x(1)) / Y
];
end
function dx = biochem(~, x, d)
% x is the state
% d is the dilution tare

% Some constants
mu_max = 0.4; % Specific grouth rate
km = 0.12; % Monod saturation constant
k1 = 0.4545; % inverse of inhibition constant
Y = 0.4; % Process Yield r1/r2

% Assume that the species input concentration 
% is zero, and that the substrate concentration
% at the input is
x2f = 4.0;

mu = mu_max * x(2) / (km + x(2) + k1 * x(2)^2);

dx = [
x(1) * (mu - d);
(x2f - x(2)) * d - (mu * x(1)) / Y
];
end