% Script for the application of the pqEDMD to the TcLab
% get the data from the file
% clear variables
load("tclab_hour.mat")
%%
dyndata = table2array(data);
tot_sam = height(dyndata);
%%
% Divide the samples, single trajectory for training last part for testing
s_ii = [1, 1350];
s_if = [1351, 1800];
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f, [2 3]), ...
	'u',dyndata(i:f, [4 5]))],s_ii,s_if);

% Training and testing
ts = 2;
tr = 1;
%%
svdEDMD = pqEDMDm(p=[2 3 4], ...
	q=[0.5 0.7 0.9 1.1 1.5 2], ...
	observable = @(p,q,l)laguerreObservable(p,q,l,0), ...
	dyn_dcp = @svdDecomposition);

dcps = svdEDMD.fit(samples(tr));
%
err = arrayfun(@(dcp)dcp.error(samples(ts)),dcps)
% calculate the best trajectories
[minerr, best] = min(err);

% Extract the best approximation
dcp = dcps(best);
% Make the prediction of the test set
svd_tst = dcp.pred_from_test(samples(ts));

%% plot it
% Define some colors
yellow  = "#FCB711";
orange  = "#F37021";
blue    = "#0089D0";
green   = "#0DB14B";
violet  = "#B711FC";
purple  = "#7021F3";
cyan			 = "#11FCB7";
springreen = "#21F370";
svd_fig = figure(1);
clf
hold on
% Plot training temperatures
plot(samples(tr).t, samples(tr).y(:,1), '.', "Color", orange, "MarkerSize", 10)
plot(samples(tr).t, samples(tr).y(:,2), '.', "Color", yellow, "MarkerSize", 10)
% Plot testing temperatures (inverted colors)
plot(samples(ts).t, samples(ts).y(:,1), '.', "Color", violet, "MarkerSize", 10)
plot(samples(ts).t, samples(ts).y(:,2), '.', "Color", purple, "MarkerSize", 10)
% Plot the heat inputs
plot(data.t, data.Q1, '--', "Color", blue, "LineWidth", 2)
plot(data.t, data.Q2, '--', "Color", green, "LineWidth", 2)
% Plot the approximation
plot(samples(ts).t, svd_tst.y(:,1), "Color", cyan, "LineWidth", 3)
plot(samples(ts).t, svd_tst.y(:,2), "Color", springreen, "LineWidth", 3)
xlabel("t [s]","Interpreter","latex")
ylabel("($^\circ$C) \& ($\%$)", "Interpreter","latex")
legend("T1","T2","T1","T2","Q1","Q2","App1","App2",Location="northoutside",Orientation="horizontal")

set(gcf,"PaperPosition",[0,0,20,10])
saveas(svd_fig,strcat("~/Documents/Seminar/figures/", "appx.fig"))
saveas(svd_fig,strcat("~/Documents/Seminar/figures/", "appx.eps"), "epsc")
