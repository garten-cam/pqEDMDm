% Script for the application of the pqEDMD to the TcLab
% get the data from the file
clear variables
load("tclab_hour.mat")
%%
dyndata = table2array(data);
tot_sam = height(dyndata);
%%
% New approach, single trajectory for training last part for testing
s_ii = [1, 1350];
s_if = [1351, 1800];
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f, [2 3]), ...
	'u',dyndata(i:f, [4 5]))],s_ii,s_if)';
str_data = struct('t', dyndata(:,1), ...
'y', dyndata(:, [2 3]), ...
'u', dyndata(:, [4 5]));

% Training and testing
ts = 2;
tr = 1;
%%
sidEDMD = pqEDMDm(p=[1 2 3],q=[0.2 0.3 0.5 0.7 1 1.2 1.4], ...
	observable = @(p,q,l)laguerreObservable(p,q,l,1), ...
dyn_dcp = @(observable, data)sidDecomposition(1, 4, observable, data));
			% dyn_dcp = @(o,s)sidDecomposition(10,6,o,s));
%
dcps = sidEDMD.fit(samples(tr));
%
err = arrayfun(@(dcp)dcp.error(samples(ts)),dcps)
% calculate the best trajectories
[minerr, best] = min(err);
% Extract the best approximation
dcp = dcps(best);
% Make the prediction of the test set
sid_tst = dcp.pred_from_test(samples(ts));

%% Filtering
% % svddcp = svdDecomposition(laguerreObservable(4,1,2,1), sm_n(tr))
filtered = dcp.filter_lift(str_data);
%% Retrain
% Divide the data into training and testing for retraining
flt_samples = arrayfun(@(i,f)[struct(...
  'y', filtered.y(i:f,:), ...
  'u', filtered.u(i:f,:))], s_ii, s_if)';

%%
% Define some colors
trainingT1 = "#332288";
testingT2  = "#44AA99";
trainingT2 = "#AA4499";
testingT1  = "#CC6677";
heatinput1 = "#332288";
heatinput2 = "#999933";
appx_T12	 = "#000000";
%
markersize = 0.1;
linewidth  = 1.0;
%% Problem statement and training set
% First figure, the trainig and testing sets
trts_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
scatter(samples(tr).t./60, samples(tr).y(:,1), markersize, "o", "MarkerEdgeColor", trainingT1, "MarkerFaceColor", trainingT1)
scatter(samples(tr).t./60, samples(tr).y(:,2), markersize, "o", "MarkerEdgeColor", trainingT2, "MarkerFaceColor", trainingT2)
% Plot the testing trajectories
scatter(samples(ts).t./60, samples(ts).y(:,1), markersize, "o", "MarkerEdgeColor", testingT1, "MarkerFaceColor", testingT1)
scatter(samples(ts).t./60, samples(ts).y(:,2), markersize, "o", "MarkerEdgeColor", testingT2, "MarkerFaceColor", testingT2)
xline(samples(tr).t(end)./60, '--k', 'LineWidth', 0.8, 'Alpha', 0.6)
xticklabels({})
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (tr)", "$T_2$ (tr)", "$T_1$ (ts)", "$T_2$ (ts)", ...
       'Interpreter', 'latex', 'Location', 'best', ...
       'Orientation', 'horizontal')
nexttile([2 1])
hold on
grid on
grid minor
plot(data.t./60, data.Q1, '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(data.t./60, data.Q2, '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")
xline(samples(tr).t(end)./60, '--k', 'LineWidth', 0.8, 'Alpha', 0.6)
legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
saveas(trts_fig,strcat("~/Documents/ECC27/figures/", "trts.fig"))
% exportgraphics(trts_fig, "~/Documents/ECC27/figures/trts.eps", "Resolution", 600);
saveas(trts_fig,strcat("~/Documents/ECC27/figures/", "trts.eps"), "epsc")
%% Plot here the first approximation
testmkrsz = 0.5;
sid_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
% Plot the testing trajectories
scatter(samples(ts).t./60, samples(ts).y(:,1), testmkrsz, "o", "MarkerEdgeColor", testingT1, "MarkerFaceColor", testingT1)
scatter(samples(ts).t./60, samples(ts).y(:,2), testmkrsz, "o", "MarkerEdgeColor", testingT2, "MarkerFaceColor", testingT2)
% Plot the approximation
plot(samples(ts).t./60, sid_tst.y(:,1), "-.", "Color", appx_T12, "LineWidth", linewidth)
plot(samples(ts).t./60, sid_tst.y(:,2), "-.", "Color", appx_T12, "LineWidth", linewidth)
xticklabels({})
ylim([25 65])
% xlim([45 60])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (ts)", "$T_2$ (ts)", "$T_1$ (appx)", "$T_2$ (appx)", ...
       'Interpreter', 'latex', 'Location', 'best', ...
       'Orientation', 'horizontal')
nexttile([2 1])
hold on
grid on
grid minor
plot(samples(ts).t./60, samples(ts).u(:,1), '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(samples(ts).t./60, samples(ts).u(:,2), '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
% xlim([45 60])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")

legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')

saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_bad.fig"))
saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_bad.eps"), "epsc")

%% Plot the filtering of all the data
filt_line = 0.5;
flt_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
scatter(samples(tr).t./60, samples(tr).y(:,1), markersize, "o", "MarkerEdgeColor", trainingT1, "MarkerFaceColor", trainingT1)
scatter(samples(tr).t./60, samples(tr).y(:,2), markersize, "o", "MarkerEdgeColor", trainingT2, "MarkerFaceColor", trainingT2)
% Plot the testing trajectories in the same color as the training
scatter(samples(ts).t./60, samples(ts).y(:,1), markersize, "o", "MarkerEdgeColor", trainingT1, "MarkerFaceColor", trainingT1)
scatter(samples(ts).t./60, samples(ts).y(:,2), markersize, "o", "MarkerEdgeColor", trainingT2, "MarkerFaceColor", trainingT2)
% Plot the filtered trajectories for training and testing
plot(data.t./60, filtered.y(:,1), "-", "Color", testingT1, "LineWidth", filt_line)
plot(data.t./60, filtered.y(:,2), "-", "Color", testingT2, "LineWidth", filt_line)
% xline(samples(tr).t(end)./60, '--k', 'LineWidth', 0.8, 'Alpha', 0.6)
xticklabels({})
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (data)", "$T_2$ (data)","","", "$T_1$ (filt)", "$T_2$ (filt)", ...
       'Interpreter', 'latex', 'Location', 'best', ...
       'Orientation', 'horizontal')
nexttile([2 1])
hold on
grid on
grid minor
plot(data.t./60, data.Q1, '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(data.t./60, data.Q2, '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")
legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
saveas(flt_fig, strcat("~/Documents/ECC27/figures/", "filt.fig"))
% exportgraphics(trts_fig, "~/Documents/ECC27/figures/trts.eps", "Resolution", 600);
saveas(flt_fig, strcat("~/Documents/ECC27/figures/", "filt.eps"), "epsc")
%% Recalculate with filtered data
% flt_sidEDMD = pqEDMDm(p=[1 2 3 4],q=[0.2 0.3 0.5 0.7 1 1.2 1.4], ...
%   observable = @(p,q,l)laguerreObservable(p,q,l,1), ...
%   dyn_dcp = @(observable, data)sidDecomposition(4, 1, observable, data));
flt_sidEDMD = pqEDMDm(p=[4 5],q=[0.2 0.3 0.5 0.7 1 1.2 1.4 2 2.3], ...
  observable = @(p,q,l)hermiteObservable(p,q,l), ...
  dyn_dcp = @(observable, data)sidDecomposition(4, 4, observable, data));
flt_dcps = flt_sidEDMD.fit(flt_samples(tr));
flt_err = arrayfun(@(flt_dcp)flt_dcp.error(flt_samples(ts)), flt_dcps)
[~, flt_idx] = min(flt_err);
flt_dcp = flt_dcps(flt_idx);
% predict from the test set
flt_tst = flt_dcp.pred_from_test(flt_samples(ts));
% Plot second appx
sid_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
% Plot the testing trajectories
scatter(samples(ts).t./60, samples(ts).y(:,1), testmkrsz, "o", "MarkerEdgeColor", testingT1, "MarkerFaceColor", testingT1)
scatter(samples(ts).t./60, samples(ts).y(:,2), testmkrsz, "o", "MarkerEdgeColor", testingT2, "MarkerFaceColor", testingT2)
% Plot the approximation
plot(samples(ts).t./60, flt_tst.y(:,1), "-.", "Color", appx_T12, "LineWidth", linewidth)
plot(samples(ts).t./60, flt_tst.y(:,2), "-.", "Color", appx_T12, "LineWidth", linewidth)
xticklabels({})
ylim([25 65])
xlim([45 60])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (ts)", "$T_2$ (ts)", "$T_1$ (appx)", "$T_2$ (appx)", ...
  'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
nexttile([2 1])
hold on
grid on
grid minor
plot(samples(ts).t./60, samples(ts).u(:,1), '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(samples(ts).t./60, samples(ts).u(:,2), '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
xlim([45 60])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")

legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')

saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_flt.fig"))
saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_flt.eps"), "epsc")