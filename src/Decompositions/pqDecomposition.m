classdef pqDecomposition
	%PQDECOMPOSITION parent class for all the decompositions.
	%
	% Returns a set of matrices characterizing a system based on the pqEDMD
	% algorithm and the provided data.
	properties % input properties
		obs % input property of every decomposition. The observable
	end
	properties % calculated properties
		A		% Evolution Matrix psi(x) = A*psi(x) + B*u
		B		% Input matrix              this one   ^
		C		% Output matrix psi(y) = C*psi(x) + D*u
		Cob	% Observables output matrix y = Cob*psi(y)
		D % Input Output matrix
		l % Number of outputs
		m % Number of inputs
		n % Number of states, in this case, number of observables
		num_obs % Number of observables
	end
	methods
		function obj = pqDecomposition(observable, data)
			%PQDECOMPOSITION Class constructor. Recieves an observable and a set
			%of trajectories/experimets/samples and  returns the A, B, and C
			%matrices of an EDMD.
			%
			% According to the observable function $\Psi$,
			%
			if nargin > 0
				
				obj.obs = observable; % saves the observable object
				
				[obs_pst, obs_fut] = obj.y_snapshots(data);
				
				if isfield(data,'u')
					obj.m = size(data(1).u, 2);
					[u_pst, u_fut] = obj.u_snapshots(data);
					obs_pst = [obs_pst u_pst];
					obs_fut = [obs_fut u_fut];
				else
					obj.m = 0;
				end
				obj.num_obs = size(obj.obs.polynomials_order, 2) + 1;
				
				% Get the regression from the snapshots
				U = obj.regression(obs_fut, obs_pst); % This is the transpose of U
				
				obj.l = obj.obs.l;
				obj.n = size(obj.obs.polynomials_order,2) + 1;
				
				% for the A, B and C properties, I have to call a method in
				% the constructor to assign them. otherwise it always
				% executes the whole get method.
				obj.A = obj.matrix_A(U);
				obj.B = obj.matrix_B(U);
				obj.C = eye(obj.num_obs); % In these decompositions psi(y)=x
				obj.Cob = obj.matrix_C;
				obj.D = zeros(obj.num_obs,max(1,obj.m));
			end
		end
		function u = regression(obj, obs_fut, obs_pst)
			% Need this simple method here to inherit later
			g = obs_pst'*obs_pst;
			a = obs_pst'*obs_fut;
			u = g\a;
		end
		function a = matrix_A(obj, U)
			% get the A matrix from the data x = Ax
			a = U(1:end-obj.m, 1:end-obj.m)';
		end
		function b = matrix_B(obj, U)
			if ~obj.m
				b = zeros(obj.n,1);
			else
				b = U(end-obj.m+1:end,1:end-obj.m)';
			end
		end
		function c = matrix_C(obj)
			orders = obj.obs.poly_base;
			order_one_obs = orders(:, 1:obj.l );
			z = sym('z',[1 obj.l],'real');
			sol = solve(order_one_obs' == z');
			if isstruct(sol)
				sol = struct2array(sol);
			end
			[AA, bb] = equationsToMatrix(sol,z);
			c = double([-bb, AA, zeros(obj.obs.l, ...
				size(obj.obs.poly_base,2) - obj.l)]);
		end
		function pred = predict(obj, y0, n_points, u)
			if nargin < 4
				u = arrayfun(@(x) {zeros(x,1)}, n_points);
			end
			% preallocate
			pred = arrayfun(@(x) struct('y', zeros(x, obj.obs.l)), n_points);
			% assign the initial condition

			obsf = obj.obs.obs_function;
			for orb = 1 : size(y0,1) % For all the initial conditions
				% Just save the initial output
				pred(orb).y(1,:) = y0(orb,:);
				for step = 2 : n_points(orb)
					% Lift the previous output
					x_prev = [1 obsf(pred(orb).y(step-1,:))];
					% Evolve the "state"
					x_post = obj.A*x_prev' + obj.B*u{orb}(step-1,:)';
					% Apply the output function
					pred(orb).y(step,:) = (obj.Cob*(obj.C*x_post + obj.D*u{orb}(step,:)'))';
				end
			end
		end
		function pred = predict_lift(obj, y0, n_points, u)
			if nargin < 4
				u = arrayfun(@(x) {zeros(x,1)}, n_points);
			end
			% preallocate
			pred = arrayfun(@(x) struct('y', zeros(x, obj.obs.l)), n_points);
			% assign the initial condition

			obsf = obj.obs.obs_function;
			for orb = 1 : size(y0,1) % For all the initial conditions
				% Just save the initial output
				pred(orb).y(1,:) = y0(orb,:);
				% linfting the initial condition
				x = [1 obsf(pred(orb).y(1,:))];
				for step = 2 : n_points(orb)
					% Lift the previous output
					% x_prev = [1 obsf(pred(orb).y(step-1,:))];
					% Evolve the "state"
					x = (obj.A*x' + obj.B*u{orb}(step-1,:)')';
					% Apply the output function
					pred(orb).y(step,:) = (obj.Cob*(obj.C*x' + obj.D*u{orb}(step,:)'))';
				end
			end
		end
		function err = error(obj, data)
			% error returns the error over all the provided testing
			% samples
			% Get the prediction to calculate error
			pred = obj.pred_from_test(data);
			% e = (y-hat(y))/y , normalized by the number of points, num of outs
			% and num of testing trajectories
			err = sum(arrayfun(@(y,hy) ... y in the training set and hat(y) in the prediction
				sum(abs(y.y-hy.y)./(abs(y.y)+eps),"all")/length(y.y)/obj.obs.l,...
				data,pred))/numel(data);
		end
		function err = abs_error(obj, data)
			% error returns the error over all the provided testing
			% samples
			% Get the prediction to calculate error
			pred = obj.pred_from_test(data);
			% e = (y-hat(y)) , normalized by the number of points, num of outs
			% and num of testing trajectories
			err = sum(arrayfun(@(y,hy) ... y in the training set and hat(y) in the prediction
				sum(abs(y.y-hy.y),"all")/length(y.y)/obj.obs.l,...
				data,pred))/numel(data);
		end
		function pred = pred_from_test(obj, data)
			% predict all the xts for this decomposition:
			% extract the initial conditions
			y0 = cell2mat(arrayfun(@(x) {x.y(1, :)}, data));
			% extract the numper of points per orbit
			np = arrayfun(@(x) size(x.y,1), data);
			% preallocate
			% I have to deal with the inputs. This was wrong....
			if isfield(data, "u")
				pred = obj.predict(y0, np, {data.u});
			else
				pred = obj.predict(y0, np);
			end
		end
		function pred = pred_from_test_lift(obj, data)
			% predict all the xts for this decomposition:
			% extract the initial conditions
			y0 = cell2mat(arrayfun(@(x) {x.y(1, :)}, data));
			% extract the numper of points per orbit
			np = arrayfun(@(x) size(x.y,1), data);
			% preallocate
			% I have to deal with the inputs. This was wrong....
			if isfield(data, "u")
				pred = obj.predict_lift(y0, np, {data.u});
			else
				pred = obj.predict_lift(y0, np);
			end
		end
		function [y_obs_pst, y_obs_fut] = y_snapshots(obj, data)
			obsrv = obj.obs.obs_function;
			y_obs_pst = cell2mat(arrayfun(@(trj) ...
				{[ones(size(trj.y,1)-2,1) obsrv(trj.y(1:end-2,:))]}, data));
			y_obs_fut = cell2mat(arrayfun(@(trj) ...
				{[ones(size(trj.y,1)-2,1) obsrv(trj.y(2:end-1,:))]}, data));
		end
		function spectrum(obj)
			figure()
			e = eig(obj.A);
			scatter(real(e),imag(e))
			hold on
			plot(sin(0:pi/100:2*pi),cos(0:pi/100:2*pi))
		end
	end
	methods (Static)
		function [u_pst, u_fut] = u_snapshots(data)
			u_pst = cell2mat(arrayfun(@(trj) ...
				{trj.u(1:end-2,:)}, data));
			u_fut = cell2mat(arrayfun(@(trj) ...
				{trj.u(2:end-1,:)}, data));
		end
	end
end
