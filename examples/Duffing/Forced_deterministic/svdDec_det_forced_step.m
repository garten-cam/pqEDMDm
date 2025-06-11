% Script for testing the Duffing equation, forced and deterministic with
% the decomposition based on the singular value decomosition with
% effectie rank
clear variables
rng(1)
% call the script that simulates the experiments

num_ics = 9;
exp = duff_det_forced_step( ...
	num_ics, ... # of initial conditions
	round(normrnd(20,2,num_ics,1)), ... # final time of all experiments
	50); % point multiplier sets the sampling time number of points = mult*tfin + 1

% The same dilemma, to normalize or not to normalize.
% exp = normalize_data(duff_exp,[-1, 2]);
% exp = duff_exp;

% Test the orthogonal pqEDMD
% ts = [2 4 5 6]; % index of training trajectories
% tr = [1 3];
tr = 1:4;
ts = 5:9;

tas_pq = pqEDMDm(p=[2 3 4 5 6], ... 
	q=[0.5 1 1.5 2], ...
	observable = @legendreObservable, ...
	dyn_dcp = @svdDecomposition); % ''  to use the ordinary least squares
dcps = tas_pq.fit(exp(tr));

%%
% Get the best performing decomposition
err = zeros(numel(dcps),1);
for decp = 1 : numel(dcps)
	err(decp) = dcps(decp).abs_error(exp(ts));
end

[er_bst, best] = min(err);
dcp = dcps(best);

appx = dcp.pred_from_test_lift(exp(ts));

tas_f = figure(1);
clf
hold on
% plot the training set
trp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b',LineWidth=2),exp(tr));
% Plot the testing set
tsp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'r',LineWidth=2),exp(ts));
% Plot the sid approximation
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=1.5),appx);
legend([trp(1), tsp(1), sia(1)],{"training", "testing", "approx"})
title({"svdDecomposition" + ...
  " p="+num2str(dcp.obs.p) + ...
  " q="+num2str(dcp.obs.q), ...
  " d=" + num2str(dcp.num_obs) + ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")
% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))