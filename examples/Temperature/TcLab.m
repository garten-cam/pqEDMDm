% Script for the application of the pqEDMD to the TcLab
% get the data from the file
dataTcL
dyndata = table2array(dyndata);
% Retrieve the samples
s_ii = linspace(1,501,6);
s_if = linspace(100,600,6);
% % some random numbers to make the trajectories of different lenghts
shift = [-8 -4 1 10 10];%%
s_ii(2:6) = s_ii(2:6) + shift;
s_if(1:5) = s_if(1:5) + shift;
%
% New approach, single trajectory for training last part for testing
% s_ii = [1, 501];
% s_if = [500, 600];
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f,4:5), ...
	'u',dyndata(i:f,2:3))],s_ii,s_if)';


% Normalize the data
sm_n = normalize_data(samples,[-1,1]);
% sm_n = samples;   

% Trainin and testing
ts = [4];
% ts = 3;
tr = [1, 2, 3, 5, 6];
% tr = 1;

sidEDMD = pqEDMDm(p=[2 3 4],q=[0.5 1 1.5], ...
	observable = @hermiteObservable, ...
	dyn_dcp = @svdDecomposition);%@(o,s)sidDecomposition(8,2,o,s));
%%
sid_dcps = sidEDMD.fit(sm_n(tr));
%
err = arrayfun(@(dcp)dcp.error(sm_n(ts)),sid_dcps);
% calculate the best trajectories
[~, best_dcp] = min(err);

% Extract the best approximation
sid_dcp = sid_dcps(best_dcp);
% Make the prediction of the test set
sid_tst = sid_dcp.pred_from_test(sm_n(ts));
% plot it
figure(1)
clf
hold on
trp_vv = arrayfun(@(sam){plot(sam.t,sam.y,'b','LineWidth',2)},sm_n(tr));
tsp_vv = arrayfun(@(sam){plot(sam.t,sam.y,'r','LineWidth',2)},sm_n(ts));
app_v1 = arrayfun(@(tst,prd){plot(tst.t,prd.y(:,1),'-.k','LineWidth',2)},sm_n(ts),sid_tst);
app_v2 = arrayfun(@(tst,prd){plot(tst.t,prd.y(:,2),'-.k','LineWidth',2)},sm_n(ts),sid_tst);
in_ht1 = arrayfun(@(sam){plot(sam.t,sam.u(:,1),'g','LineWidth',2)},sm_n);
in_ht2 = arrayfun(@(sam){plot(sam.t,sam.u(:,2),'m','LineWidth',1.5)},sm_n);
