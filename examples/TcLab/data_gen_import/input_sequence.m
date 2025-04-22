function q = input_sequence(duration, delta_t, step_interval)
	%INPUT_SEQUENCE Generates a sequence of inputs for the TcLab
	%
	% duration: duration in seconds of the experiment
	% delta_t: sampling time, or time step in seconds
	% step_interval: number of delta_t that a step signal stays activated
	
	np = ceil(duration/delta_t); % Total number of points

	q = zeros(np,1);
	for segment = 1 : step_interval : np
		q(segment:min(segment+step_interval-1,np)) = randi([0 90]);
	end
end

