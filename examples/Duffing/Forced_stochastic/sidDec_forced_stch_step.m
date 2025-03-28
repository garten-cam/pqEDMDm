% Script for testing the Duffing equation, forced and deterministic with
% the decomposition based on the singular value decomosition with
% effectie rank
clear variables
% call the script that simulates the experiments

exp_forced_stch_step

% The same dilemma, to normalize or not to normalize.
% exp = normalize_data(duff_exp,[-1, 1]);
exp = duff_exp;
% Test the orthogonal pqEDMD
% ts = [2 4 5 6]; % index of training trajectories
% tr = [1 3];
tr = 1:7;
ts = 8:9;

tas_pq = pqEDMDm(p=[2 3 4 5], ... [2 3 4 5], ...
	q=[0.5 1 1.5 2 2.5], ...
	observable = @legendreObservable, ...
	dyn_dcp = @(sy,ob)sidDecomposition(5,1,sy,ob)); % ''  to use the ordinary least squares
dcps = tas_pq.fit(exp(tr));
%%
% Get the best performing decomposition from the error of the testing set
err = arrayfun(@(dcpi)dcpi.error(exp(ts)),dcps);

[er_bst, best] = min(err);
dcp = dcps(best);

appx = dcp.pred_from_test(exp(ts));

tas_f = figure(1);
clf
hold on
% plot the deterministic samples
det = arrayfun(@(ex)plot(ex.y_det(:,1),ex.y_det(:,2),'g',LineWidth=3),exp);
% plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b',LineWidth=2),exp(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r',LineWidth=2),exp(ts));
% Plot the sid approximation
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=1.5),appx);
legend([det(1), trp(1), tsp(1), sia(1)],{"determ","training", "testing", "approx"})
title({"sidDecomposition" + ...
  " p="+num2str(dcp.obs.p) + ...
  " q="+num2str(dcp.obs.q), ...
  " d=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n), ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")
% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))