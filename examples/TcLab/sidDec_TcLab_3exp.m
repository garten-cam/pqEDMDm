% New attempt with the old data. Why did I feel I had to use another set, I
% have no idea.
clear variables
% Get the 4 available datasets
data{1} = importdata("tclab_hourE1.mat");
data{2} = importdata("tclab_hourE2.mat");
data{3} = importdata("tclab_hourE3.mat");
%%
resample = 1;
samples = cellfun(@(dti)[struct(...
  't',dti.t(1:resample:end,1), ...
  'y',[dti.T1(1:resample:end) dti.T2(1:resample:end)], ...
  'u',[dti.Q1(1:resample:end) dti.Q2(1:resample:end)])],data)';
% Define the training and testing sets
tr = [1 2];
ts = 3;
% Define the qp sweep over all the decompositions
sidEDMD = pqEDMDm(p=[1 2 3 4],q=[0.25 0.5 1 1.25], ...
	observable = @(p,q,l)laguerreObservable(p,q,l,1), ...
	dyn_dcp = @(o,d)sidDecomposition(2, 4, o, d));
% Fit the thing
dcps = sidEDMD.fit(samples(tr));
%
err = arrayfun(@(dcp)dcp.error(samples(ts)),dcps)
%
[minerr, best] = min(err);
%
% Extract the best approximation
dcp = dcps(best);
% Make the prediction of the test set
sid_tst = dcp.pred_from_test(samples(ts));
%% Filter
filtered = dcp.filter_lift(samples);
% Check the spectral radious
flt_err = sum(abs(samples(ts).y-filtered(ts).y)./(abs(samples(ts).y)+eps),"all")/length(filtered(ts).y)/dcp.obs.l;
% make a small test from zero initial conditions
zero_ic = samples;
for sam = 1:length(samples)
  zero_ic(sam).y(1,:) = [0 0];
end
zero_flt = dcp.filter_lift(zero_ic);

%% Plot here the first approximation
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
testmkrsz = 0.5;
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
scatter(samples(tr(1)).t./60, samples(tr(1)).y(:,1), markersize, "o", "MarkerEdgeColor", trainingT1, "MarkerFaceColor", trainingT1)
scatter(samples(tr(1)).t./60, samples(tr(1)).y(:,2), markersize, "o", "MarkerEdgeColor", trainingT2, "MarkerFaceColor", trainingT2)
xticklabels({})
ylim([15 70])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$", "$T_2$", ...
  'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
nexttile([2 1])
hold on
grid on
grid minor
plot(samples(tr(1)).t./60, samples(tr(1)).u(:,1), '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(samples(tr(1)).t./60, samples(tr(1)).u(:,2), '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time [min]", Interpreter="latex")
legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
% saveas(trts_fig,strcat("~/Documents/ECC27/figures/", "trts.fig"))
% % exportgraphics(trts_fig, "~/Documents/ECC27/figures/trts.eps", "Resolution", 600);
% saveas(trts_fig,strcat("~/Documents/ECC27/figures/", "trts.eps"), "epsc")

%% Plot the approximation on the test set
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
ylim([15 80])
% xlim([45 60])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (ts)", "$T_2$ (ts)", "$T_1$ (appx)", "$T_2$ (appx)", ...
  'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')

title({"$p$="+num2str(dcp.obs.p) + ...
  " $q$="+num2str(dcp.obs.q) + ...
  " $d$=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n) + ...
  " $\epsilon$="+num2str(minerr,2)}, Interpreter='latex')

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

% saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_bad.fig"))
% saveas(sid_fig,strcat("~/Documents/ECC27/figures/", "sid_appx_bad.eps"), "epsc")

%% Plot the filtered test set
filt_line = 0.5;
flt_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
scatter(samples(ts).t./60, samples(ts).y(:,1), testmkrsz, "o", "MarkerEdgeColor", testingT1, "MarkerFaceColor", testingT1)
scatter(samples(ts).t./60, samples(ts).y(:,2), testmkrsz, "o", "MarkerEdgeColor", testingT2, "MarkerFaceColor", testingT2)
% Plot the filtered trajectories for training and testing
plot(samples(ts).t./60, filtered(ts).y(:,1), "-", "Color", trainingT1, "LineWidth", filt_line)
plot(samples(ts).t./60, filtered(ts).y(:,2), "-", "Color", trainingT2, "LineWidth", filt_line)

xticklabels({})
ylim([15 80])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (ts)", "$T_2$ (ts)","$T_1$ (filt)", "$T_2$ (filt)", ...
  'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
title({"$p$="+num2str(dcp.obs.p) + ...
  " $q$="+num2str(dcp.obs.q) + ...
  " $d$=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n) + ...
  " $\epsilon$="+num2str(flt_err,2)}, Interpreter='latex')
nexttile([2 1])
hold on
grid on
grid minor
plot(samples(ts).t./60, samples(ts).u(:,1), '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(samples(ts).t./60, samples(ts).u(:,2), '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")
legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')

saveas(flt_fig, strcat("~/Documents/ECC27/figures/", "filt.fig"))
saveas(flt_fig, strcat("~/Documents/ECC27/figures/", "filt.eps"), "epsc")
%%
% Test the warmup
wrm_stp = 60;
wrm_plt = 6;
wrmup_ex = dcp.filter_warm_up(samples(ts), wrm_stp);
wrmup_err = sum(abs(samples(ts).y(1:wrm_stp*wrm_plt,:)-wrmup_ex.y(1:wrm_stp*wrm_plt,:))./(abs(samples(ts).y(1:wrm_stp*wrm_plt,:))+eps),"all")/length(samples(ts).y(1:wrm_stp*wrm_plt,:))/dcp.obs.l;
 
filt_line = 0.5;
wrmup_fig = figure('Units', 'centimeters', 'Position', [0, 0, 18, 8]);
clf
tiledlayout("vertical","TileSpacing","none","Padding","tight");
% Plot the training trajectories
nexttile([3 1])
hold on
grid on
grid minor
scatter(samples(ts).t(1:wrm_stp*wrm_plt)./60, samples(ts).y(1:wrm_stp*wrm_plt,1), testmkrsz, "o", "MarkerEdgeColor", testingT1, "MarkerFaceColor", testingT1)
scatter(samples(ts).t(1:wrm_stp*wrm_plt)./60, samples(ts).y(1:wrm_stp*wrm_plt,2), testmkrsz, "o", "MarkerEdgeColor", testingT2, "MarkerFaceColor", testingT2)
% Plot the filtered trajectories for training and testing
plot(samples(ts).t(1:wrm_stp*wrm_plt)./60, wrmup_ex.y(1:wrm_stp*wrm_plt,1), "-", "Color", trainingT1, "LineWidth", filt_line)
plot(samples(ts).t(1:wrm_stp*wrm_plt)./60, wrmup_ex.y(1:wrm_stp*wrm_plt,2), "-", "Color", trainingT2, "LineWidth", filt_line)

xticklabels({})
% ylim([15 80])
ylabel("$T$ [$^\circ$C]", "Interpreter", "latex")
legend("$T_1$ (ts)", "$T_2$ (ts)","$T_1$ (filt)", "$T_2$ (filt)", ...
  'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
title({"$p$="+num2str(dcp.obs.p) + ...
  " $q$="+num2str(dcp.obs.q) + ...
  " $d$=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n) + ...
  " $\epsilon$="+num2str(wrmup_err,2)}, Interpreter='latex')
nexttile([2 1])
hold on
grid on
grid minor
plot(samples(ts).t(1:wrm_stp*wrm_plt)./60, samples(ts).u(1:wrm_stp*wrm_plt,1), '-.', "Color", heatinput1, "LineWidth", linewidth)
plot(samples(ts).t(1:wrm_stp*wrm_plt)./60, samples(ts).u(1:wrm_stp*wrm_plt,2), '-.', "Color", heatinput2, "LineWidth", linewidth)
ylim([0 95])
ylabel("Power ($\%$)", "Interpreter", "latex")
yticks([0 50 90])
xlabel("Time $t$ [min]", Interpreter="latex")
legend("$Q_1$", "$Q_2$", 'Interpreter', 'latex', 'Location', 'best', ...
  'Orientation', 'horizontal')
saveas(wrmup_fig, strcat("~/Documents/ECC27/figures/", "wrmup.fig"))
saveas(wrmup_fig, strcat("~/Documents/ECC27/figures/", "wrmup.eps"), "epsc")