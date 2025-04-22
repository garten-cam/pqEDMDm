% Script to run the EDMD tutorial using a TcLab.

clear variables
% Set plotting colors
yellow = '#FCB711';
orange = '#F37021';
blue   = '#0089D0';
green  = '#0DB14B';
rng(1) % for reproducibility

% Heat inputs parameters
duration   = 3600;
delta_t		 = 2; % log every 2 seconds
step_dt_Q1 = 120; % Change every 2 minutes
step_dt_Q2 = 180; % Change every minure and a half
% Generate the two input sequences
Q1_sec = input_sequence(duration, delta_t, step_dt_Q1);
Q2_sec = input_sequence(duration, delta_t, step_dt_Q2); % Change every 20 seconds

%
np = ceil(duration/delta_t); % number of data points

% Preallocate the values
data = table( ...
	'Size',[np,5], ...
	'VariableTypes',{'double','double','double','double','double'}, ...
	'VariableNames',{'t','T1','T2','Q1','Q2'});

tclab
figure(1)
for ts = 0 : np - 1
	tic
	%
	% Get the hater values from the sequence
	Q1 = Q1_sec(ts+1);
	Q2 = Q2_sec(ts+1);
	% Set the heaters
	h1(Q1)
	h2(Q2)
	% Read the temperatures
	T1 = T1C();
	T2 = T2C();
	% Log the data
	data.t(ts+1)  = ts*delta_t;
	data.T1(ts+1) = T1;
	data.T2(ts+1) = T2;
	data.Q1(ts+1) = Q1;
	data.Q2(ts+1) = Q2;

	% Plot the temperatures
	clf
	hold on
	plot(data.t(1:ts+1), data.T1(1:ts+1), '.', ...
		'Color', hex2rgb(orange), 'MarkerSize', 10);
	plot(data.t(1:ts+1), data.T2(1:ts+1), '.', ...
		'Color', hex2rgb(yellow), 'MarkerSize', 10);
	% Plot the heat inputs
	plot(data.t(1:ts+1), data.Q1(1:ts+1), '--', ...
		'Color', hex2rgb(blue), 'LineWidth', 2);
	plot(data.t(1:ts+1), data.Q2(1:ts+1), '--', ...
		'Color', hex2rgb(green), 'LineWidth', 2);
	ylabel("Temperature [$^\circ$C] and Heat input [$\%$]", 'Interpreter','latex')
	xlabel("t [s]","Interpreter","latex")
	legend("T1","T2","Q1","Q2")
	drawnow

	% Wait until the next time to sample
	pause(delta_t - mod(toc,delta_t))
end
% end the experiment
h1(0)
h2(0)
