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
tfin = 30;
n_points = 20*tfin;
%
% Create the input signal. for simplicity, it will be a step response,
% i.e., a constant input in the range [0.1 - 1]
in = (0.28-0.08)*rand(num_ics, 1)+0.08;


% select now a training and a testing set
tr = [1 2 4 6 7 8 9 10];
ts = [3 5]; % The second argumet is the num of test cycles

%%
% test wit noise.
mu_max = .4*ones(num_ics,1);% + .4/10 * rand(num_ics, 1);
% 
exp_stch = repmat(struct('y', zeros(n_points + 1, 2), ...
  'u', zeros(n_points + 1, 1),... forcing signals
  't', zeros(n_points + 1, 1)), num_ics,1);


odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
  [exp_stch(orb).t, exp_stch(orb).y] = ode23s(@(t,x)biochem_stch(t,x, ...
      in(orb),...
      mu_max(orb)),...
    0:tfin/n_points:tfin, ...
    ics(orb,:), ...
    odeSettings);
  % exp_stch(orb).y = exp_stch(orb).y + exp_stch(orb).y/10 .* rand(size(exp_stch(orb).y));
  exp_stch(orb).y = exp_stch(orb).y + normrnd(0,.01,size(exp_stch(orb).y));
  exp_stch(orb).u = in(orb)*ones(size(exp_stch(orb).t));
  % cycles(orb).u = cos(in(orb)*cycles(orb).t);
end
%%
% Plot the thing according to time
figure(1)
time_resp = tiledlayout(2,5,"TileSpacing","tight");
for smp = 1:numel(exp_stch)
  nexttile(smp)
  hold on
  plot(exp_stch(smp).t, exp_stch(smp).y(:,1), 'b')
  plot(exp_stch(smp).t, exp_stch(smp).y(:,2), 'r')
  if smp == 1
    legend("$x_1$ Bio","$x_2$ Sub","Interpreter","latex", "location", "best")
  end
end
% arrayfun(@(smp)plot(smp.t,smp.y(:,1),'b'),exp_stch)
% arrayfun(@(smp)plot(smp.t,smp.y(:,2),'r'),exp_stch)
%%
% Robust decomposition
pen_pq_sidR = pqEDMDm( ...
    p = [2 3],...
    q = [0.5 1],...
    observable = @legendreObservable,...
    dyn_dcp = @sidDecomposition...
    );

% Calculate the decomposition
pen_stchR = pen_pq_sidR.fit(exp_stch(tr));
%%
% Test the decompositions
err_sidR = zeros(numel(pen_stchR),1);
for decp = 1 : numel(pen_stchR)
  err_sidR(decp) = pen_stchR(decp).error(exp_stch(ts));
end
% where is the min?
[~, best_sidR] = min(err_sidR)
%
st_sp_sidR = figure(2);
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

% Predict from the test 
pen_p_stchR = pen_stchR(best_sidR).pred_from_test(exp_stch(ts));

% Plot the result
figure(2)
for cyc = 1 : numel(ts)
    tsa = plot(pen_p_stchR(cyc).y(:,1), pen_p_stchR(cyc).y(:,2),'-.k',LineWidth=2);
end
legend([trp, tsp, tsa],{'training','testing','EDMD'})
title("sid Decomposiotion")
%%
pen_pq_sidO = pqEDMDm( ...
    p = [2 3 4],...
    q = [0.5 1],...
    observable = @legendreObservable,...
    dyn_dcp = @sidOlsDecomposition);
 
% Calculate the decomposition
pen_stchO = pen_pq_sidO.fit(exp_stch(tr));
%%
% Test the decompositions
err_sidO = zeros(numel(pen_stchO),1);
for decp = 1 : numel(pen_stchO)
  err_sidO(decp) = pen_stchO(decp).abs_error(exp_stch(ts));
end
% where is the min?
[minerr_sid,best_sidO] = min(err_sidO)
%%
st_sp_sidO = figure(3);
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

% Predict from the test 
pen_p_stchO = pen_stchO(best_sidO).pred_from_test(exp_stch(ts));

% Plot the result
figure(3)
for cyc = 1 : numel(ts)
    tsa = plot(pen_p_stchO(cyc).y(:,1), pen_p_stchO(cyc).y(:,2),'-.k',LineWidth=2);
end
legend([trp, tsp, tsa],{'training','testing','sidEDMD'})
title("sidO Decomposiotion")


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