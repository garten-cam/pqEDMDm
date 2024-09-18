% Script for the application of the pqEDMD to the TcLab
% get the data from the file
clear variables
dataTcL
dyndata = table2array(dyndata);
tot_sam = height(dyndata);
% Retrieve the samples
segments = 5;
s_ii = floor(1:tot_sam/segments:tot_sam);
% %s_ii = floor(linspace(1,600,7));
s_if = [s_ii(2:end)-1 tot_sam];
% % some random numbers to make the trajectories of different lenghts
% shift = [-8 -4 1 10 10];%%
% s_ii(2:6) = s_ii(2:6) + shift;
% s_if(1:5) = s_if(1:5) + shift;
%
% New approach, single trajectory for training last part for testing
% s_ii = [1, 501];
% s_if = [500, 600];
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f,4), ...
	'u',dyndata(i:f,[2 5]))],s_ii,s_if)';


% Normalize the data
% sm_n = normalize_data(samples,[-1,1]);
sm_n = samples;   

% Trainin and testing
tr = 1:segments;
allcomb = nchoosek(tr,2);
ts = allcomb(randi(length(allcomb),1),:);
tr(ts) = [];
% ts = 3;
% tr = [1, 4, 5, 6, 7, 8];
% tr = 1;
%%
sidEDMD = pqEDMDm(p=[4 5 6 7 8],q=[0.5 1 1.5 2], ...
	observable = @hermiteObservable, ...
	dyn_dcp = @(o,s)sidDecomposition(4,3,o,s));
		% dyn_dcp = @svdDecomposition);
%
dcps = sidEDMD.fit(sm_n(tr));
%
err = arrayfun(@(dcp)dcp.abs_error(sm_n(ts)),dcps)
% calculate the best trajectories
[minerr, best] = min(err);

% Extract the best approximation
dcp = dcps(best);
% Make the prediction of the test set
sid_tst = dcp.pred_from_test(sm_n(ts));
% plot it
figure(1)
clf
hold on
trp_vv = arrayfun(@(sam){plot(sam.t,sam.y,'b','LineWidth',2)},sm_n(tr));
tsp_vv = arrayfun(@(sam){plot(sam.t,sam.y,'r','LineWidth',2)},sm_n(ts));
app_v1 = arrayfun(@(tst,prd){plot(tst.t,prd.y(:,1),'-.k','LineWidth',2)},sm_n(ts),sid_tst);
%
% app_v2 = arrayfun(@(tst,prd){plot(tst.t,prd.y(:,2),'-.k','LineWidth',2)},sm_n(ts),sid_tst);
in_ht1 = arrayfun(@(sam){plot(sam.t,sam.u(:,1),'g','LineWidth',2)},sm_n);
in_ht2 = arrayfun(@(sam){plot(sam.t,sam.u(:,2),'m','LineWidth',1.5)},sm_n);

