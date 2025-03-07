% Script for testing the Duffing equation, forced and deterministic with
% the decomposition based on the singular value decomosition with
% effectie rank
clear variables
rng(34)
% call the script that simulates the experiments
num_ics = 9;
exp = duff_det_forced_step( ...
	num_ics, ... # of initial conditions
	round(normrnd(20,2,num_ics,1)), ... # final time of all experiments
	5); % point multiplier sets the sampling time number of points = mult*tfin + 1

% The same dilemma, to normalize or not to normalize.
% exp = normalize_data(duff_exp,[-1, 1]);
% exp = duff_exp;
% Test the orthogonal pqEDMD
% ts = [2 4 5 6]; % index of training trajectories
% tr = [1 3];
tr = [1, 3];
ts = [2, 4:9];

tas_pq = pqEDMDm(p=5, ... [2 3 4 5], ... 
	q=1.5, ... [0.5 1 1.5 2]
	observable = @legendreObservable, ...
	dyn_dcp = @(obs,data)sidDecomposition(50,1,obs,data)); % ''  to use the ordinary least squares
dcps = tas_pq.fit(exp(tr));
%%
% Get the best performing decomposition from the error of the testing set
err = arrayfun(@(dcpi)dcpi.abs_error(exp(ts)),dcps);

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
sia = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'-.k',LineWidth=1.5),appx([1,2,3,4,5,6,7]));
legend([trp(1), tsp(1), sia(1)],{"training", "testing", "approx"})
title({"sidDecomposition" + ...
  " p="+num2str(dcp.obs.p) + ...
  " q="+num2str(dcp.obs.q), ...
  " d=" + num2str(dcp.num_obs) + ...
  " $n$="+num2str(dcp.n), ...
  " $\epsilon$="+num2str(er_bst)},Interpreter="latex")
% axis([-1 1 -1 1])
%saveas(tas_f,strcat(figpath, "tr_ts.png"))
%saveas(tas_f,strcat(figpath, "approx.png"))
%%
tas_t = figure(2);
clf
tile_t = tiledlayout("vertical","TileSpacing","tight","Padding","tight");
nexttile
hold on
arrayfun(@(ex){plot(ex.t,ex.y,'b',LineWidth=2)},exp(tr));
nexttile
hold on
arrayfun(@(ex){plot(ex.t,ex.y,'r',LineWidth=1.5)},exp(ts([1,2])));
arrayfun(@(exd,exa){plot(exd.t,exa.y,'-.k',LineWidth=1.5)},exp(ts([1,2])),appx([1,2]))
