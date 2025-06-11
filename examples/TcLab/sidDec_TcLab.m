% Script for the application of the pqEDMD to the TcLab
% get the data from the file
 clear variables
load("tclab_hour.mat")
%%
dyndata = table2array(data);
tot_sam = height(dyndata);
%%
% Retrieve the samples
% segments = 2;
% s_ii = floor(1:tot_sam/segments:tot_sam);
% %s_ii = floor(linspace(1,600,7));
% s_if = [s_ii(2:end)-1 tot_sam];
% % some random numbers to make the trajecto$ries of different lenghts
% shift = [-8 -4 1 10 10];%%
% s_ii(2:6) = s_ii(2:6) + shift;
% s_if(1:5) = s_if(1:5) + shift;
%
% New approach, single trajectory for training last part for testing
s_ii = [1, 1350];
s_if = [1351, 1800];
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f, [2 3]), ...
	'u',dyndata(i:f, [4 5]))],s_ii,s_if)';


% Normalize the data
% sm_n = normalize_data(samples,[0,1]);
sm_n = samples; 

% Trainin and testing
% tr = 1:segments;
% allcomb = nchoosek(tr,2);
% ts = allcomb(randi(length(allcomb),1),:);
% tr(ts) = [];
ts = 2;
% tr = [1, 4, 5, 6, 7, 8];
tr = 1;
%%
sidEDMD = pqEDMDm(p=[1 2 3 4],q=[0.5 1 2 3], ...
	observable = @(p,q,l)laguerreObservable(p,q,l,1), ...
dyn_dcp = @(data, observable)sidDecomposition(4, 4, data, observable));
			% dyn_dcp = @(o,s)sidDecomposition(10,6,o,s));
%
dcps = sidEDMD.fit(sm_n(tr));
%
err = arrayfun(@(dcp)dcp.error(sm_n(ts)),dcps)
% calculate the best trajectories
[minerr, best] = min(err);

% Extract the best approximation
dcp = dcps(best);
% Make the prediction of the test set
sid_tst = dcp.pred_from_test(sm_n(ts));
% Define some colors
yellow  = "#FCB711";
orange  = "#F37021";
blue    = "#0089D0";
green   = "#0DB14B";
violet  = "#B711FC";
purple  = "#7021F3";
cyan			 = "#11FCB7";
springreen = "#21F370";

sid_fig = figure(2);
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
plot(samples(ts).t, sid_tst.y(:,1), "Color", cyan, "LineWidth", 3)
plot(samples(ts).t, sid_tst.y(:,2), "Color", springreen, "LineWidth", 3)
xlabel("t [s]","Interpreter","latex")
ylabel("($^\circ$C) \& ($\%$)", "Interpreter","latex")
legend("T1","T2","T1","T2","Q1","Q2","App1","App2",Location="northoutside",Orientation="horizontal")

set(gcf,"PaperPosition",[0,0,20,10])
saveas(sid_fig,strcat("~/Documents/Seminar/figures/", "sid_appx_bad.fig"))
saveas(sid_fig,strcat("~/Documents/Seminar/figures/", "sid_appx_bad.eps"), "epsc")
