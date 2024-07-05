% Script for the application of the pqEDMD to the TcLab
% get the data from the file
dataTcL
dyndata = table2array(dyndata);
% Retrieve the samples
s_ii = linspace(1,501,6);
s_if = linspace(100,600,6);
% some random numbers to make the trajectories of different lenghts
shift = [-8 -4 1 10 10];%%
s_ii(2:6) = s_ii(2:6) + shift;
s_if(1:5) = s_if(1:5) + shift;
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
	'y',dyndata(i:f,4:5), ...
	'u',dyndata(i:f,2:3))],s_ii,s_if)';

% Normalize the data
sm_n = normalize_data(samples,[-1,1]);

% Trainin and testing
ts = [3, 5];
tr = [1, 2, 4, 6];

sidEDMD = pqEDMDm(p=[4 5],q=[1.2 1.5], ...
	observable = @legendreObservable, ...
	dyn_dcp = @svdDecomposition);
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
in_set = arrayfun(@(sam){plot(sam.t,sam.u,'g','LineWidth',2)},sm_n)


%%
% plots....
figure(1)
clf
hold on
for cyc = tr
	trp = plot(samples(cyc).t,samples(cyc).y,'b');
	plot(samples(cyc).t,samples(cyc).u,'g')
end
for cyc = 1:length(ts)
	tsp = plot(samples(ts(cyc)).t,samples(ts(cyc)).y,'r');
	tst = plot(samples(ts(cyc)).t,testR(cyc).y,'k');
	plot(samples(ts(cyc)).t,samples(ts(cyc)).u,'g')
end
title('Robust')
figure(2)
clf
hold on
for cyc = tr
	trp = plot(samples(cyc).t,samples(cyc).y,'b');
	plot(samples(cyc).t,samples(cyc).u,'g')
end
for cyc = 1:length(ts)
	tsp = plot(samples(ts(cyc)).t,samples(ts(cyc)).y,'r');
	tst = plot(samples(ts(cyc)).t,testO(cyc).y,'k');
	plot(samples(ts(cyc)).t,samples(ts(cyc)).u,'g')
end
title('Ols')


