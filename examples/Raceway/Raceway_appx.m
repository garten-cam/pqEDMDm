% Load the data
load("raceway_data.mat")

% Identification

% Define the training and the testing sets
tr = [1 2 3 4]; % some for trainig
ts = [5 6]; % some for testing
% Define the pqEDMD wrapper
rw_edmd = pqEDMDm(p=[2 3 4],q=[0.5 1 1.5 2], ...
	observable = @(a,b,c)laguerreObservable(a,b,c,1), ...
dyn_dcp = @svdDecomposition);
%
% Fit all the p-q pairs
dcps = rw_edmd.fit(exps(tr));
%
% Calculate the error of all the p-q pairs
err = arrayfun(@(dcp)dcp.error(exps(ts)),dcps);
% Which is the best?
[minerr, best] = min(err);
%
dcp = dcps(best);
% Calculate the prediction of the test set from the best pq combo
tst = dcp.pred_from_test(exps(ts));
%
% Plot the outputs
trs_fig = figure(1);
clf
trs_tile = tiledlayout(7,1,"TileSpacing","tight","Padding","tight");
% title(trs_tile, 'First Experiment',Interpreter='latex')
nexttile
hold on
exp_cb_plt_1 = plot(exps(1).t, exps(1).y(:,2));
ylabel('Cb [g/m\textsuperscript{3}]',Interpreter='latex')
xsecondarylabel(Visible='off')
% legend('training',"Interpreter","latex",Location="southwest")
title('(a) biomass concentration',Interpreter='latex')
grid on
set(gca,'XTicklabel',[])
%
nexttile
exp_do_plt_1 = plot(exps(1).t, exps(1).y(:,3));
ylabel('DO [g/m\textsuperscript{3}]',Interpreter='latex')
hold on
yldo = yline(DOrs(1),'r');
legend(yldo,'$O^\star$',Interpreter='latex')
xsecondarylabel(Visible='off')
grid on
set(gca,'XTicklabel',[])
title('(b) dissolved oxygen',Interpreter='latex')
%
nexttile
exp_ph_plt_1 = plot(exps(1).t, exps(1).y(:,4));
ylabel('pH [-]',Interpreter='latex')
hold on
ylph = yline(pHrs(1),'r');
legend(ylph,'$pH^\star$',Interpreter='latex')
xsecondarylabel(Visible='off')
grid on
set(gca,'XTicklabel',[])
title('(c) pH',Interpreter='latex')
%
nexttile
exp_tmp_plt_1 = plot(exps(1).t,exps(1).u(:,1));
ylabel('T [$^\circ$C]',Interpreter='latex')
xsecondarylabel(Visible='off')
grid on
set(gca,'XTicklabel',[])
title('(d) water temperature',Interpreter='latex')
%
nexttile
exp_rad_plt_1 = plot(exps(1).t,exps(1).u(:,2));
ylabel('rad [W/m\textsuperscript{2}]',Interpreter='latex')
xsecondarylabel(Visible='off')
grid on
set(gca,'XTicklabel',[])
title('(e) solar radiation',Interpreter='latex')
%
nexttile
exp_air_plt_1 = plot(exps(1).t,exps(1).u(:,5));
ylabel('Air [m\textsuperscript{3}/s]',Interpreter='latex')
xsecondarylabel(Visible='off')
grid on
set(gca,'XTicklabel',[])
title('(f) air injection',Interpreter='latex')
%
nexttile
exp_co2_plt_1 = plot(exps(1).t,exps(1).u(:,6));
ylabel('CO$_2$ [m\textsuperscript{3}/s]',Interpreter='latex')
xlabel('Time [h]',Interpreter='latex')
grid on
xsecondarylabel(Visible='off')
title('(g) CO$_2$ injection',Interpreter='latex')
%
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% Plot the results. Just the testing set
tst_fig = figure(2);
clf
stt_tile = tiledlayout(3,2,"TileSpacing","tight","Padding","tight");
% title(stt_tile,'Testing Set',Interpreter='latex')
%
nexttile
hold on
exp_cb_plt_5 = plot(exps(5).t, exps(5).y(:,2));
tst_cb_plt_5 = plot(exps(5).t, tst(1).y(:,2),'-.','LineWidth',1);
ylabel('Cb [g/m\textsuperscript{3}]',Interpreter='latex')
xsecondarylabel(Visible='off')
legend('test','appx',"Interpreter","latex",Location="best")
title('(a) biomass day 5',Interpreter='latex')
grid on
set(gca,'XTicklabel',[])
%
nexttile
hold on
set(gca, 'ColorOrderIndex',1)
exp_cb_plt_6 = plot(exps(5).t, exps(5).y(:,2));
tst_cb_plt_6 = plot(exps(5).t, tst(1).y(:,2),'-.','LineWidth',1);
xsecondarylabel(Visible='off')
title('(b) biomass day 6',Interpreter='latex')
grid on
set(gca,'XTicklabel',[])
set(gca,'YTickLabel',[])
%
nexttile
hold on
set(gca, 'ColorOrderIndex',1)
exp_do_plt_5 = plot(exps(5).t, exps(5).y(:,3));
tst_do_plt_5 = plot(exps(5).t, tst(1).y(:,3),'-.','LineWidth',1);
ylabel('DO [g/m\textsuperscript{3}]',Interpreter='latex')
xsecondarylabel(Visible='off')
title('(c) oxygen day 5',Interpreter='latex')
grid on
set(gca,'XTicklabel',[])
%
nexttile
hold on
set(gca, 'ColorOrderIndex',1)
exp_do_plt_6 = plot(exps(6).t, exps(6).y(:,3));
tst_do_plt_6 = plot(exps(6).t, tst(2).y(:,3),'-.','LineWidth',1);
xsecondarylabel(Visible='off')
title('(d) oxygen day 6',Interpreter='latex')
grid on
set(gca,'XTicklabel',[])
set(gca,'YTickLabel',[])
%
nexttile
hold on
set(gca, "ColorOrderIndex",1)
exp_ph_plt_5 = plot(exps(5).t, exps(5).y(:,4));
tst_ph_plt_5 = plot(exps(5).t, tst(1).y(:,4),'-.','LineWidth',1);
ylabel('pH [-]',Interpreter='latex')
xlabel('Time [h]', Interpreter='latex')
xsecondarylabel(Visible='off')
title('(e) pH day 5',Interpreter='latex')
grid on
%
nexttile
hold on
set(gca, "ColorOrderIndex",1)
exp_ph_plt_6 = plot(exps(6).t, exps(6).y(:,4));
tst_ph_plt_6 = plot(exps(6).t, tst(2).y(:,4),'-.','LineWidth',1);
title('(f) pH day 6',Interpreter='latex')
grid on
% set(gca,'XTicklabel',[])
set(gca,'YTickLabel',[])
xlabel('Time [h]',Interpreter='latex')
xsecondarylabel(Visible='off')
