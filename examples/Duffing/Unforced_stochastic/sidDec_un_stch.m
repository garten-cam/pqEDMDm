% Script for testing the Duffing equation, unforced, and deterministic with
% the basic pqDecomposition
clear variables
% call the script that generates the experiments

exp_un_stch_meas;

% To normalize or not. In this case is not neccesary, leaving it here for
% reference
% exp = normalize_data(duff_exp,[-1,1]);

% "useless" copy, basicatly just renaming in case we want to test the
% normalized case. 
exp = duff_exp;

% Test the orthogonal pqEDMD
ts = [2 4 5 6]; % index of training trajectories
tr = [1 3];

tas_pq = pqEDMDm(p=[3 4 5], ... p=[2 3 4 5]
	q=[2 2.5], ...
	observable = @legendreObservable, ...
	dyn_dcp = @(ob,sy)sidOlsDecomposition(1,1,ob,sy)); % '' to use the ordinary least squares
dcps = tas_pq.fit(exp(tr));
%%
% Get the best performing decomposition
err = zeros(numel(dcps),1);
for decp = 1 : numel(dcps)
	err(decp) = dcps(decp).abs_error(exp(ts));
end

[er_bst, best] = min(err);
dcp = dcps(best);

appx = dcp.pred_from_test(exp(ts));

tas_f = figure(1);
clf
hold on
% pot the deterministic samples
det = arrayfun(@(ex)plot(ex.y_det(:,1),ex.y_det(:,2),'g',LineWidth=3),exp);
% plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b',LineWidth=2),exp(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r',LineWidth=2),exp(ts));
% Plot the sid approximation
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=1.5),appx);
legend([det(1), trp(1), tsp(1), sia(1)],{"determ","training", "testing", "approx"})
title({"sidDecomposition" + ...
  " $p$="+num2str(dcp.obs.p) + ...
  " $q$="+num2str(dcp.obs.q), ...
  " $d$=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n), ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))