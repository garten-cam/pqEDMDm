% Script for testing the Duffing equation, forced and deterministic with
% the decomposition based on the singular value decomosition with
% effectie rank
clear variables
rng(1234)
% call the script that simulates the experiments
num_ics = 6;
exp = duff_stch_forced_cos( ...
	num_ics, ... # of initial conditions
	round(normrnd(15,2,num_ics,1)), ... # final time of all experiments
	0.1, ... Standard deviation for the noise addition
	20); % point multiplier sets the sampling time number of points = mult*tfin + 1

% Test the orthogonal pqEDMD
% ts = [2 4 5 6]; % index of training trajectories
% tr = [1 3];
tr = 1:num_ics-2;
ts = num_ics-1:num_ics;
%%
tas_pq = pqEDMDm(p=[2 3 4], ... [2 3 4 5], ...
	q=[0.5 1 1.5 2], ...
	observable = @legendreObservable, ...
	dyn_dcp = @(sy,ob)sidDecomposition(2,4,sy,ob)); % ''  to use the ordinary least squares
dcps = tas_pq.fit(exp(tr));
%
% Get the best performing decomposition from the error of the testing set
err = arrayfun(@(dcpi)dcpi.abs_error(exp(ts)),dcps)

[er_bst, best] = min(err);
dcp = dcps(best);

appx = dcp.pred_from_test(exp(ts));
%
tas_f = figure(1);
tl = tiledlayout('flow','TileSpacing','tight','Padding','tight');
% Plot the deterministic samples
for ex = 1 : numel(exp)
	nexttile(ex);
	hold on
	det = plot(exp(ex).t, exp(ex).y_det, 'g', LineWidth=2);
end
% Plot the training set
for tri = tr
	nexttile(tri)
	trp = plot(exp(tri).t, exp(tri).y, 'b', LineWidth=1);
end
% Plot the testing set
for tsi = ts
	nexttile(tsi)
	tsp = plot(exp(tsi).t, exp(tsi).y, 'r', LineWidth=1);
end
% Plot the approximation
for apxi = 1 : numel(appx)
	nexttile(ts(apxi))
	apxp = plot(exp(ts(apxi)).t, appx(apxi).y, '-.k');
end

legend([det(1), trp(1), tsp(1), apxp(1)],{"determ","training", "testing", "approx"})
title(tl,{"sidDecomposition" + ...
  " p="+num2str(dcp.obs.p) + ...
  " q="+num2str(dcp.obs.q), ...
  " d=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n), ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")
% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))