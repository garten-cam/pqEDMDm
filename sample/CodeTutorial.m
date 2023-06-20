% Script for producing the examples of the tutorial paper.
% clear variables
rng(1)
% Start with the basics. The duffing equation in all of its
% representations
%%
%----------------------------------------------------------------------
% Dataset for the simulations
%----------------------------------------------------------------------
%  Set the same initial conditions for the three cases
num_ics = 6; % Number of initial conditions for the test
ics_width = 4; % ics range width
ics = ics_width*rand(num_ics,2) - ics_width/2;

%%

% Define the simulation structure for all the cases. The same initial
% conditions.
t_fin = 20;
n_points = 10*t_fin + 1;
sim_param = struct('ICs',ics,...
    'uP',zeros(2,num_ics),...
    'nP',[0 0],...
    'TFin',t_fin,...
    'NPoints',n_points);
no_stop = false;
%%
% Duffing parameters for undamped, hard stifness and soft stiffness
dup.alpha = 1;
dup.beta = 1;
dup.delta = 0;
% Calling the Orbits class without ODE, input function, or noise function requires a name value pair for the simulation structure
undamped_samples = Orbits(@(t,x,n,u)DuffEqODE(t,x,dup,n,u),'Solver','ode23s','simparam',sim_param,'stop',no_stop);
%  GASo
dup.alpha = 1;
dup.beta = 1;
dup.delta = 0.5;

% Calling the Orbits class without ODE, input function, or noise function requires a name value pair for the simulation structure
gas_samples = Orbits(@(t,x,n,u)DuffEqODE(t,x,dup,n,u),'Solver','ode23s','simparam',sim_param,'stop',no_stop);
%%
% Two AS
dup.alpha = -1;
dup.beta = 1;
dup.delta = 0.5;

% Calling the Orbits class without ODE, input function, or noise function requires a name value pair for the simulation structure
as_samples = Orbits(@(t,x,n,u)DuffEqODE(t,x,dup,n,u),'Solver','ode23s','simparam',sim_param,'stop',no_stop);

%%
% Figures
num_to_plot = 2;
axdisp = 3;
lw=1.1;
lz = 12; %labelsize
samples_to_plot = randperm(num_ics,num_to_plot);
%%
phase_lines_und = gobjects(num_ics,1);
duffplot = figure(1);
clf
subplot(2,3,1)
hold on
for orbit = 1 :num_ics
    phase_lines_und(orbit) = plot(undamped_samples.clean_orbits(orbit).SV(:,1),undamped_samples.clean_orbits(orbit).SV(:,2),'b','LineWidth',lw);
end
% g = scatter(0,0,20,'g','filled');
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
% legend([h(1),g],{'$\alpha=1$, $\beta=1$, $\delta=0$','Stable'},'Location',"northwest",'Interpreter',"latex")
legend(phase_lines_und(1),'$\alpha=1$, $\beta=1$, $\delta=0$','Location',"northwest",'Interpreter',"latex")
axis([-axdisp axdisp -axdisp axdisp])

% Plot against time. Not all of them because it gets crouded, just
% one.

time_lines = gobjects(num_to_plot,2);
subplot(2,3,4)
hold on

for orbit = 1 : num_to_plot
    time_lines(orbit,1) = plot(undamped_samples.clean_orbits(samples_to_plot(orbit)).Time,undamped_samples.clean_orbits(samples_to_plot(orbit)).SV(:,1),'-m');
    time_lines(orbit,2) = plot(undamped_samples.clean_orbits(samples_to_plot(orbit)).Time,undamped_samples.clean_orbits(samples_to_plot(orbit)).SV(:,2),'-g');
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x$','Interpreter',"latex","FontSize",lz)
legend([time_lines(1,1),time_lines(1,2)],{'$x_1$','$x_2$'},'Location',"northwest",'Interpreter',"latex")

%%
phase_lines_gas = gobjects(num_ics,1);
figure(1);
subplot(2,3,2)
hold on
for orbit = 1 : num_ics
    phase_lines_gas(orbit) = plot(gas_samples.clean_orbits(orbit).SV(:,1),gas_samples.clean_orbits(orbit).SV(:,2),'b','LineWidth',lw);
end
% g = scatter(0,0,20,'g','filled');
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend(phase_lines_gas(1),'$\alpha=1$, $\beta=1$, $\delta=0.5$','Location',"northwest",'Interpreter',"latex")
axis([-axdisp axdisp -axdisp axdisp])

% Plot against time. Not all of them because it gets crouded, just
% one.
time_lines_gas = gobjects(num_to_plot,2);
subplot(2,3,5)
hold on
for orbit = 1 : num_to_plot
    time_lines_gas(orbit,1) = plot(gas_samples.clean_orbits(samples_to_plot(orbit)).Time,gas_samples.clean_orbits(samples_to_plot(orbit)).SV(:,1),'-m');
    time_lines_gas(orbit,2) = plot(gas_samples.clean_orbits(samples_to_plot(orbit)).Time,gas_samples.clean_orbits(samples_to_plot(orbit)).SV(:,2),'-g');
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x$','Interpreter',"latex","FontSize",lz)
legend([time_lines_gas(1,1),time_lines_gas(1,2)],{'$x_1$','$x_2$'},'Location',"northeast",'Interpreter',"latex")

%%
phase_lines_as = gobjects(num_ics,1);
figure(1);
subplot(2,3,3)
hold on
for orbit = 1 : num_ics
    phase_lines_as(orbit) = plot(as_samples.clean_orbits(orbit).SV(:,1),as_samples.clean_orbits(orbit).SV(:,2),'b','LineWidth',lw);
end
% g = scatter(0,0,20,'g','filled');
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend(phase_lines_as(1),'$\alpha=-1$, $\beta=1$, $\delta=0.5$','Location',"northwest",'Interpreter',"latex")
axis([-axdisp axdisp -axdisp axdisp])

% Plot against time. Not all of them because it gets crouded, just
% one.
time_lines_as = gobjects(num_to_plot,2);
subplot(2,3,6)
hold on
for orbit = 1 : num_to_plot
    time_lines_as(orbit,1) = plot(as_samples.clean_orbits(samples_to_plot(orbit)).Time,as_samples.clean_orbits(samples_to_plot(orbit)).SV(:,1),'-m');
    time_lines_as(orbit,2) = plot(as_samples.clean_orbits(samples_to_plot(orbit)).Time,as_samples.clean_orbits(samples_to_plot(orbit)).SV(:,2),'-g');
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x$','Interpreter',"latex","FontSize",lz)
legend([time_lines_as(1,1),time_lines_as(1,2)],{'$x_1$','$x_2$'},'Location',"northeast",'Interpreter',"latex")
set(gcf,'Paperposition',[0 0 15 10])

%  
%%
% Start from the simple case GAS. I need to compare the  
% Get the EDMD for the nostop
[~,far_ic] = max(vecnorm(ics,1,2));
% [~,near_ic] = min(vecnorm(ics,1,2));
ts_index = 1:1:num_ics;
tr_index = ts_index([far_ic]);
ts_index(tr_index) = [];
gas_nostop_approx_full = pqEDMD(gas_samples.clean_orbits, [3 4 6], [0.7],...
    'polynomial','Legendre','polyParam',[2 3],'method','OLS',...
    'tr_ts',struct('tr_index',tr_index,'ts_index',ts_index),'normalization',false);
gas_nostop_approx = gas_nostop_approx_full.best_pqEDMD
gas_nostop_test_approx = struct('Time', cell(length(ts_index),1),...
                                'SV', cell(length(ts_index),1));
evol_function = matlabFunction(gas_nostop_approx.h_psi(gas_nostop_approx.VVFO.polynomial_base*gas_nostop_approx.U),'var',{sym("x",[1 2])});
for orbit = 1 : length(ts_index) % Outer loop, for all the available test orbits
    gas_nostop_test_approx(orbit).Time = gas_samples.clean_orbits(ts_index(orbit)).Time;
    for datapoint  = 1 : length(gas_samples.clean_orbits(ts_index(orbit)).Time) % Inner loop over all the datapoints
        if datapoint == 1
            gas_nostop_test_approx(orbit).SV(datapoint,:) = gas_samples.clean_orbits(ts_index(orbit)).SV(datapoint,:); % If it is the first iteration, save the real initial condition
        else
%             txn_approx(orbit).LagMult
            gas_nostop_test_approx(orbit).SV(datapoint,:) = evol_function(normalize(gas_nostop_test_approx(orbit).SV(datapoint-1,:),'center',gas_nostop_approx.center,'scale',gas_nostop_approx.scale)).*gas_nostop_approx.scale+gas_nostop_approx.center;
        end
    end
end

% empty objects for the legends
%%
test_gas = figure(2);
clf
subplot(2,2,[1 3])
hold on
for training_orb = 1 : length(tr_index)
    gas_training_plot = plot(gas_samples.clean_orbits(tr_index(training_orb)).SV(:,1),gas_samples.clean_orbits(tr_index(training_orb)).SV(:,2),'r','LineWidth',1.5*lw);
end
for test_orb = 1 : length(ts_index)
    gas_testing_plot = plot(gas_samples.clean_orbits(ts_index(test_orb)).SV(:,1),gas_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',lw);
    gas_testing_app_plot = plot(gas_nostop_test_approx(test_orb).SV(:,1),gas_nostop_test_approx(test_orb).SV(:,2),'-.k','LineWidth',lw);
end
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([gas_training_plot,gas_testing_plot,gas_testing_app_plot],...
 {'tr','ts','approx'},'Location',"northwest",'Interpreter',"latex")
for test_orb = 1 : length(ts_index)
    subplot(2,2,2)
    hold on
    smoothPlot = plot(gas_samples.clean_orbits(ts_index(test_orb)).Time,gas_samples.clean_orbits(ts_index(test_orb)).SV(:,1),'b','LineWidth',lw);
%     noisePlot = scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,1),scatter_size,scatter_shape(test_orb),'g');
    approPlot = plot(gas_nostop_test_approx(test_orb).Time,gas_nostop_test_approx(test_orb).SV(:,1),'-.k','LineWidth',lw);
    subplot(2,2,4)
    hold on
    smoothPlotx2 = plot(gas_samples.clean_orbits(ts_index(test_orb)).Time,gas_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',lw);
%     scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,2),scatter_size,scatter_shape(test_orb),'g')
    approPlotx2 = plot(gas_nostop_test_approx(test_orb).Time,gas_nostop_test_approx(test_orb).SV(:,2),'-.k','LineWidth',lw);
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([smoothPlot,approPlot],{'$x_1$ ts','$x_1$ approx'},'Location',"northeast",'Interpreter',"latex")
subplot(2,2,2)
ylabel('$x_1$','Interpreter',"latex","FontSize",lz)
legend([smoothPlotx2,approPlotx2],{'$x_2$ ts','$x_2$ approx'},'Location',"northeast",'Interpreter',"latex")
set(gcf,'PaperPosition',[0 0 20 10])


%%
% create the structures for the pqEDMD for the undamped case
% for the undamped, take the extremes as the trainiing set
[~,far_ic] = max(vecnorm(ics,1,2));
[~,near_ic] = min(vecnorm(ics,1,2));
ts_index = 1:1:num_ics;
tr_index = ts_index([far_ic near_ic]);
ts_index(tr_index) = [];
% ts_index(tr_index) = [];
% Test the current pqEDMD 
undamped_nostop_approx_full = pqEDMD(undamped_samples.clean_orbits, [3 4], [1.1 inf],...
    'polynomial','Hermite','polyParam',[2 3],'method','OLS',...
    'tr_ts',struct('tr_index',tr_index,'ts_index',ts_index),'normalization',false);

undamped_nostop_approx = undamped_nostop_approx_full.best_pqEDMD

% plot the results of the test set
undamped_test_approx = struct('Time', cell(length(ts_index),1),...
                                'SV', cell(length(ts_index),1));
evol_function = matlabFunction(undamped_nostop_approx.h_psi(undamped_nostop_approx.VVFO.polynomial_base*undamped_nostop_approx.U),'var',{sym("x",[1 2])});
for orbit = 1 : length(ts_index) % Outer loop, for all the available test orbits
    undamped_test_approx(orbit).Time = undamped_samples.clean_orbits(ts_index(orbit)).Time;
    for datapoint  = 1 : length(undamped_samples.clean_orbits(ts_index(orbit)).Time) % Inner loop over all the datapoints
        if datapoint == 1
            undamped_test_approx(orbit).SV(datapoint,:) = undamped_samples.clean_orbits(ts_index(orbit)).SV(datapoint,:); % If it is the first iteration, save the real initial condition
        else
%             txn_approx(orbit).LagMult
            undamped_test_approx(orbit).SV(datapoint,:) = evol_function(normalize(undamped_test_approx(orbit).SV(datapoint-1,:),'center',undamped_nostop_approx.center,'scale',undamped_nostop_approx.scale)).*undamped_nostop_approx.scale+undamped_nostop_approx.center;
        end
    end
end

% empty objects for the legends

test_undamped = figure(4);
clf
for training_orb = 1 : length(tr_index)
     subplot(2,2,[1 3])
    hold on
    undamped_training_plot = plot(undamped_samples.clean_orbits(tr_index(training_orb)).SV(:,1),undamped_samples.clean_orbits(tr_index(training_orb)).SV(:,2),'r','LineWidth',1.5*lw);
end
for test_orb = 1 : length(ts_index)
    subplot(2,2,[1 3])
    hold on
    undamped_testing_plot = plot(undamped_samples.clean_orbits(ts_index(test_orb)).SV(:,1),undamped_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',1.5*lw);
    undamped_testing_app_plot = plot(undamped_test_approx(test_orb).SV(:,1),undamped_test_approx(test_orb).SV(:,2),'-.k','LineWidth',lw);
end
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([undamped_training_plot,undamped_testing_plot,undamped_testing_app_plot],...
 {'tr','ts','approx'},'Location',"northwest",'Interpreter',"latex")
for test_orb = 1 : length(ts_index)
    subplot(2,2,2)
    hold on
    smoothPlot = plot(undamped_samples.clean_orbits(ts_index(test_orb)).Time,undamped_samples.clean_orbits(ts_index(test_orb)).SV(:,1),'b','LineWidth',lw);
%     noisePlot = scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,1),scatter_size,scatter_shape(test_orb),'g');
    approPlot = plot(undamped_test_approx(test_orb).Time,undamped_test_approx(test_orb).SV(:,1),'-.k');
    subplot(2,2,4)
    hold on
    smoothPlotx2 = plot(undamped_samples.clean_orbits(ts_index(test_orb)).Time,undamped_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',lw);
%     scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,2),scatter_size,scatter_shape(test_orb),'g')
    approPlot = plot(undamped_test_approx(test_orb).Time,undamped_test_approx(test_orb).SV(:,2),'-.k');
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([smoothPlot,approPlot],{'$x_1$ ts','$x_1$ approx'},'Location',"northeast",'Interpreter',"latex")
subplot(2,2,2)
ylabel('$x_1$','Interpreter',"latex","FontSize",lz)
legend([smoothPlotx2,approPlotx2],{'$x_2$ ts','$x_2$ approx'},'Location',"northeast",'Interpreter',"latex")
set(gcf,'PaperPosition',[0 0 20 10])

%%

ts_index = 1:1:num_ics;
tr_index = ts_index([1 3]);
ts_index(tr_index) = [];
as_nostop_approx_full = pqEDMD(as_samples.clean_orbits, [5], [1.0],...
    'polynomial','Legendre','polyParam',[2 3],'method','ML',...
    'tr_ts',struct('tr_index',tr_index,'ts_index',ts_index),'normalization',false);
as_nostop_approx = as_nostop_approx_full.best_pqEDMD
as_nostop_test_approx = struct('Time', cell(length(ts_index),1),...
                                'SV', cell(length(ts_index),1));
evol_function = matlabFunction(as_nostop_approx.h_psi(as_nostop_approx.VVFO.polynomial_base*as_nostop_approx.U),'var',{sym("x",[1 2])});
for orbit = 1 : length(ts_index) % Outer loop, for all the available test orbits
    as_nostop_test_approx(orbit).Time = as_samples.clean_orbits(ts_index(orbit)).Time;
    for datapoint  = 1 : length(as_samples.clean_orbits(ts_index(orbit)).Time) % Inner loop over all the datapoints
        if datapoint == 1
            as_nostop_test_approx(orbit).SV(datapoint,:) = as_samples.clean_orbits(ts_index(orbit)).SV(datapoint,:); % If it is the first iteration, save the real initial condition
        else
%             txn_approx(orbit).LagMult
            as_nostop_test_approx(orbit).SV(datapoint,:) = evol_function(normalize(as_nostop_test_approx(orbit).SV(datapoint-1,:),'center',as_nostop_approx.center,'scale',as_nostop_approx.scale)).*as_nostop_approx.scale+as_nostop_approx.center;
        end
    end
end

% empty objects for the legends

test_as = figure(5);
clf
subplot(2,2,[1 3])
hold on
for training_orb = 1 : length(tr_index)
    as_training_plot = plot(as_samples.clean_orbits(tr_index(training_orb)).SV(:,1),as_samples.clean_orbits(tr_index(training_orb)).SV(:,2),'r','LineWidth',1.5*lw);
end
for test_orb = 1 : length(ts_index)
    as_testing_plot = plot(as_samples.clean_orbits(ts_index(test_orb)).SV(:,1),as_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',lw);
    as_testing_app_plot = plot(as_nostop_test_approx(test_orb).SV(:,1),as_nostop_test_approx(test_orb).SV(:,2),'-.k','LineWidth',lw);
end
xlabel('$x_1$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([as_training_plot,as_testing_plot,as_testing_app_plot],...
 {'tr','ts','approx'},'Location',"northwest",'Interpreter',"latex")
for test_orb = 1 : length(ts_index)
    subplot(2,2,2)
    hold on
    smoothPlot = plot(as_samples.clean_orbits(ts_index(test_orb)).Time,as_samples.clean_orbits(ts_index(test_orb)).SV(:,1),'b','LineWidth',lw);
%     noisePlot = scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,1),scatter_size,scatter_shape(test_orb),'g');
    approPlot = plot(as_nostop_test_approx(test_orb).Time,as_nostop_test_approx(test_orb).SV(:,1),'-.k','LineWidth',lw);
    subplot(2,2,4)
    hold on
    smoothPlotx2 = plot(as_samples.clean_orbits(ts_index(test_orb)).Time,as_samples.clean_orbits(ts_index(test_orb)).SV(:,2),'b','LineWidth',lw);
%     scatter(system(ts_index(test_orb)).Time,system(ts_index(test_orb)).SV(:,2),scatter_size,scatter_shape(test_orb),'g')
    approPlotx2 = plot(as_nostop_test_approx(test_orb).Time,as_nostop_test_approx(test_orb).SV(:,2),'-.k','LineWidth',lw);
end
xlabel('$k$','Interpreter',"latex","FontSize",lz)
ylabel('$x_2$','Interpreter',"latex","FontSize",lz)
legend([smoothPlot,approPlot],{'$x_1$ ts','$x_1$ approx'},'Location',"northeast",'Interpreter',"latex")
subplot(2,2,2)
ylabel('$x_1$','Interpreter',"latex","FontSize",lz)
legend([smoothPlotx2,approPlotx2],{'$x_2$ ts','$x_2$ approx'},'Location',"northeast",'Interpreter',"latex")
set(gcf,'PaperPosition',[0 0 20 10])

%%
