% Script for testing the Duffing equation, unforced, and deterministic with
% the svdDecomposition

% call the script that generates the experiments

exp_unforced_det;

% To normalize or not. In this case is not neccesary, leaving it here for
% reference
% exp = normalize_data(duff_exp,[-1,1]);

% "useless" copy, basically just renaming in case we want to test the
% normalized case. In fact, the normalized case is worse in this system
exp = duff_exp;

% Test the orthogonal pqEDMD
ts = [2 4 5 6]; % index of training trajectories
tr = [1 3];

tas_pq = pqEDMDm(p=[3 4 5 6], ... 
	q=[0.5 1 1.5 2], ...
	observable = @legendreObservable, ...
	dyn_dcp = @(obs,sys)sidOlsDecomposition(4,0,obs,sys));
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
% plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b',LineWidth=2),exp(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r',LineWidth=2),exp(ts));
% Plot the sid approximation
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=1.5),appx);
legend([trp(1), tsp(1), sia(1)],{"training", "testing", "approx"}, ...
  Interpreter="latex",Location="northeast")
title({"sidDecomposition" + ...
  " $p$="+num2str(dcp.obs.p) + ...
  " $q$="+num2str(dcp.obs.q), ...
  " $d$=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n), ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")
%axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))