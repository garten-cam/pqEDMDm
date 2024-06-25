% Script for the application of the pqEDMD to the TcLab
% get the data from the file
dataTcL
dyndata = table2array(dyndata);
% Retrieve the samples
s_ii = linspace(1,501,6);
s_if = linspace(100,600,6);
% some random numbers to make the trajectories of different lenghts
shift = [-8 -4 1 10 10];%%
s_ii(2:6) = s_ii(2:6) + rand_shift;
s_if(1:5) = s_if(1:5) + rand_shift;
samples = arrayfun(@(i,f)[struct('t',dyndata(i:f,1), ...
  'y',dyndata(i:f,4:5), ...
  'u',dyndata(i:f,2:3))],s_ii,s_if)';

ts = [3, 5];
tr = [1, 2, 4, 6];

sidREDMD = pqEDMDm(p=[2 3],q=[0.6 1.2], ...
  observable = @legendreObservable, ...
  dyn_dcp = @sidDecomposition);
sidOEDMD = pqEDMDm(p=[2 3],q=[0.6 1.2], ...
  observable = @legendreObservable, ...
  dyn_dcp = @sidOlsDecomposition);
%%
sidRT = sidREDMD.fit(samples(tr));
sidOT = sidOEDMD.fit(samples(tr));
%
err_R = zeros(numel(sidRT),1);
for decp = 1 : numel(sidRT)
  err_R(decp) = sidRT(decp).error(samples(ts));
end
err_O = zeros(numel(sidOT),1);
for decp = 1 : numel(sidOT)
  err_O(decp) = sidOT(decp).error(samples(ts));
end

% calculate the best trajectories 
[~, minR] = min(err_R);
[~, minO] = min(err_O);

testR = sidRT(minR).pred_from_test(samples(ts));
for r = 1 : length(testR)
  testR(r).t = samples(ts(r)).t;
end
testO = sidOT(minO).pred_from_test(samples(ts));
for r = 1 : length(testO)
  testO(r).t = samples(ts(r)).t;
end
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


