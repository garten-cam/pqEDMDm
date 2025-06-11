classdef sidDecomposition < svdDecomposition
	%SIDDECOMPOSITION subsytem identification EDMD
	properties
		K % Kalman gain of the solution
		fb % number of Hankel blocks for the Future calculations
		pb % Past blocks.
	end
	properties (Hidden)
		det
		unforced
		n_cols % Number of columns in each sample
	end
	methods
		function obj = sidDecomposition(fb, pb, observable, data)
			%SIDDECOMPOSITION computes the subspace identification of the EDMD
			if nargin > 0
				obj.obs = observable;
				[yeval_sys, u_sys] = obj.yu_eval(data);

				obj.unforced = ~isfield(data, 'u');
				obj.det = ~(pb > 0);
				if obj.unforced
					obj.m = 0;
				else
					obj.m = size(data(1).u,2);
				end
				obj.num_obs = size(obj.obs.polynomials_order, 2) + 1; % number of outputs, number of observables
				obj.fb = fb; % Hankel blocks for the future
				% we need at least one past block
				obj.pb = max(pb,1);
				% obj.pb = pb;
				obj.l = obj.obs.l; % Number of outputs
				% Hankelize
				% Ysid = [Yp;Yf]
				Ysid = cellfun(@(sys){obj.block_hankel(sys', obj.fb, obj.pb)},yeval_sys)'; % the same 2hb
				% Usid=[Up;Uf]
				if obj.unforced
					% create empty matrices for the input when the system is
					% unforced
					Usid = cellfun(@(sys){zeros(0,width(sys))},Ysid);
				else
					Usid = cellfun(@(x) ...
						{obj.block_hankel(x', obj.fb, obj.pb)},u_sys)';
				end
				% Get the system
				[obj.A, obj.B, obj.C, obj.D, obj.K, obj.n] = obj.ABCDKn(Ysid, Usid);
				% back to the state
				obj.Cob = obj.matrix_C;
			end
		end
		function [a, b, c, d, k, n] = ABCDKn(obj, Ysid, Usid)
			[U, S] = obj.compute_svd(Ysid, Usid);
			ns = obj.num_obs:sum(diag(S)>1e-9);
			% Optimice for n
			nv = arrayfun(@(n)obj.n_cost(n, U, S, Ysid, Usid), ns);
			[~, midx] = min(nv);
			% It is cheaper to test all and select the best than to optimize.
			n = ns(midx);
			% calculate everything that the "optimization" calculates inside for
			% the best of the orders n. Maybe save and index? Nah...
			[~, res, a, b, c, d] = obj.n_cost(n, U, S, Ysid, Usid);
			SigWE = res*res'/(size(res,2)-1);
			QQ = SigWE(1:obj.n,1:obj.n);
			RR = SigWE(obj.n+1:obj.n+obj.num_obs,obj.n+1:obj.n+obj.num_obs);
			SS = SigWE(1:obj.n,obj.n+1:obj.n+obj.num_obs);
			[~,k,~,~] = idare(a',c',QQ,RR,SS); % Kalman filter ARE
			k = k'; % Kalman gain
			if obj.m == 0 % redefine so they are zeroes instead of empties
				b = zeros(n,1);
				d = zeros(obj.num_obs,1);
			end
		end
		function [u, s] = compute_svd(obj, Ysid, Usid)
			% From the hankel matrices, calculate the oblique projection
			[y_f, u_f, w_p] = obj.fut_pst_mat(Ysid, Usid);
			if obj.unforced
				if obj.det
					WOW = cell2mat(y_f);
				else
					WOW = cell2mat(cellfun(@(yi,wpi){obj.zProj_x(yi,wpi)},y_f,w_p));
				end
			else
				Oi = cellfun(@(yi,ui,wi){obj.obliqueP(yi,ui,wi)},y_f, u_f, w_p);
				WOW = cell2mat(cellfun(@(oii,ui){obj.zProjOx(oii,ui)},Oi,u_f));
			end
			[u,s,~] = svd(WOW,"econ");
		end
		function [cost, res, a, b, c, d] = n_cost(obj, n, U, S, Ysid, Usid)
			gam_inv = pinv(U(:,1:n) * sqrtm(S(1:n, 1:n)));
			% future and past matrices
			[y_f, u_f, w_p] = obj.fut_pst_mat(Ysid, Usid);
			z_i = obj.get_zi(y_f, u_f, w_p);
			% Calculate the state sequences
			x = cellfun(@(zi){gam_inv*zi},z_i);
			% Solve a regression like the original EDMD
			% Get the divided matrices
			[xp, xf, yiip, ufp, uff] = obj.edmd_like_division(x, y_f, u_f);
			% [xf;uf] = [A,B][xp;up];
			ab_lhs = [cell2mat(xf);cell2mat(uff)];
			ab_rhs = [cell2mat(xp);cell2mat(ufp)];
			ab = obj.regression(ab_lhs', ab_rhs')';
			a = ab(1:n,1:n);
			if obj.m == 0
				b = zeros(n,0);
			else
				b = ab(1:n,n+1:end);
			end
			% Solve for C
			% yiip = Cxf ---> makes no sense. But finally made it work
			c_lhs = cell2mat(yiip);
			c_rhs = cell2mat(xf);
			c = obj.regression(c_lhs', c_rhs')';
			d = zeros(height(c),width(b));
			% resc = c_lhs - cc*c_rhs;
			res = [cell2mat(xf);cell2mat(yiip)] ...
				- [a b;c d]*ab_rhs;
			cost = sum(abs(res),'all');
		end % function
		function z_i = get_zi(obj, y_f, u_f, w_p)
			if obj.unforced
				if obj.det
					z_i = y_f;
				else
					z_i = cellfun(@(yi,wi){obj.zProj_x(yi, wi)},y_f, w_p); % Uf
				end
			else
				z_i = cellfun(@(yi,ui,wi){obj.obliqueP(yi,ui,wi)},y_f, u_f, w_p);
			end
		end
		function pred = predict(obj, y0, n_points, u)
			if nargin < 4
				u = arrayfun(@(x) {zeros(x,1)}, n_points);
			end
			% preallocate
			pred = arrayfun(@(x) struct('sv', zeros(x, obj.n), ... this will save the state
				'y', zeros(x, obj.obs.l)), n_points);
			obsf  = obj.obs.obs_function;
			% C inverse to get the state
			Ci = pinv(obj.C);
			% Simulate
			for orb = 1 : size(y0,1) % for all initial conditions
				pred(orb).y(1,:) = y0(orb,:); % Save the initial output
				% Save the initial state
				% pred(orb).sv(1,:) = Cinv*(obsf(pred(orb).y(1,:))'-obj.D*u{orb}(1,:)');
				pred(orb).sv(1,:) = Ci*([1; obsf(pred(orb).y(1,:))']);
				for step = 2 : n_points(orb)
					% Lift the previous output
					lft = [1 obsf(pred(orb).y(step-1,:))];
					% Lifted output to state space
					x_prev = Ci*(lft');
					% Evolve
					x_post = obj.A*x_prev + obj.B*u{orb}(step,:)';
					% Save the state
					pred(orb).sv(step,:) = x_post';
					% Save the output
					pred(orb).y(step,:) = obj.Cob*(obj.C*x_post + obj.D*u{orb}(step-1,:)');
				end
			end
		end
		function [yev, u] = yu_eval(obj, data)
			% lifting
			obsrv = obj.obs.obs_function;
			yev = arrayfun(@(sys) {[ones(size(sys.y,1),1) obsrv(sys.y)]}, data);
			% extract the inputs into a cell
			if ~isfield(data,'u')
				u = [];
			else
				u = {data.u}';
			end
		end
		function xPzAy = obliqueP(obj, X, Y, Z)
			% Oi = X/_{Y}(Z) the oblique projection of X, onto
			% the row space of Z along Y. For that,
			% Oi = X/_{Y}(Z) = (X/oY)*pinv(Z/oY)*Z,
			% where the (X/oY) notation is the projection of Y onto the
			% orthogonal complement of Y
			% X/(oY)
			xPOy = obj.zProjOx(X, Y);
			% Z/(oY) inv
			zPOy_i = pinv(obj.zProjOx(Z, Y));
			% Ready for the projection of X onto Z along Y
			xPzAy = xPOy*zPOy_i*Z;
		end
		function [y_f, u_f, w_p] = fut_pst_mat(obj, Ysid, Usid)
			y_f = cellfun(@(y_i){y_i(end-obj.fb*(obj.num_obs)+1:end,:)},Ysid);
			u_f = cellfun(@(u_i){u_i(end-obj.fb*obj.m+1:end,:)},Usid);
			w_p = cellfun(@(ui,yi){[ui(1:obj.pb*obj.m,:);yi(1:obj.pb*obj.num_obs,:)]},Usid, Ysid);
		end
		function [x_p, x_f, yii_p, u_fp, u_ff] = edmd_like_division(obj, x, y_f, u_f)
			% edmd_like_division p: past f: future
			x_p = cellfun(@(xi){xi(:,1:end-1)},x);
			x_f = cellfun(@(xi){xi(:,2:end)},x);
			yii_p = cellfun(@(yi,xi){yi(1:obj.num_obs,1:width(xi)-1)},y_f,x);
			u_fp = cellfun(@(ui,xi){ui(1:obj.m,1:width(xi)-1)},u_f,x);
			u_ff = cellfun(@(ui,xi){ui(1:obj.m,2:width(xi))},u_f,x);
		end% Function
	end
	methods (Static)
		function zpox = zProjOx(z, x)
			% zProjOx is the orthogonal projection of z onto the row pace
			% of the orthogonal complement of x.
			% zpox = z/ox
			zpox = z - z*x'*pinv(x*x')*x;
		end
		function zp_x = zProj_x(z, x)
			% zProj_x is the orthogonal projection of z onto the row pace
			% of x. z/_x
			zp_x = z*x'*pinv(x*x')*x;
		end
		function H = block_hankel(mat, fb, pb)
			% block_hankel turns a matrix into a Hankel Matrix with fb+pb number
			% of Hankel blocks.
			% mat: is the matrix to hankelize
			% fb : is the number of "future" blocks
			% pb : is the number of "past" blocks
			H = cell2mat(arrayfun(@(blk) ...
				{mat(:, blk : blk + (width(mat) - fb - pb + 1) - 1)}', ...
				1:(fb + pb))');
		end
	end % end static methods
end % end class

% function [b, d] = BD(obj, y_eval, u_sys)
% 			% Calculates the B matrix of the system. I will assume that D
% 			% is zero, and see if I finally get a system that does not
% 			% diverge.
% 			if isempty(u_sys)
% 				b = [];
% 				d = [];
% 				return
% 			end
% 			% The last approach goes back to the data, and the data, is noisy.
% 			% instead, use the more effective approach Delgado
%
% 			% max samples
% 			max_sam = max(cellfun(@(y_i)[size(y_i,1)],y_eval));
% 			ac = obj.fullGamma(max_sam);
% 			% Now, we need the kron prod of the input and the AC prod
% 			% summed over all the available inputs in a trajectory
% 			% preallocate
% 			ukrac = cellfun(@(u_i){arrayfun(@(i){kron(zeros(1,size(u_i,2)),obj.C)},1:size(u_i,1))},u_sys);
% 			% populate
% 			for u_i = 1 : length(u_sys)
% 				for u_blk = 1 : length(y_eval{u_i})
% 					krcb = zeros([size(ukrac{u_i}{u_blk}), u_blk]);
% 					for step = 1 : u_blk
% 						krcb(:,:,step) = kron(u_sys{u_i}(step,:),ac{u_blk-step+1});
% 					end
% 					ukrac{u_i}{u_blk} = sum(krcb,3);
% 				end
% 			end
% 			% Ok, I have the kroneckers, now I need the matrix
% 			% All the trajectories have the same CA;;;CA^k matrix. but,
% 			% truncated to the number of samples. so.........
% 			% I have length(y_sys) number of blocks. in the right hand side
% 			% Preallocate the blocks...
% 			%------------------
% 			rhsblk = cellfun(@(y_i){zeros(numel(y_i(2:end,:)),obj.n*(length(y_eval)+obj.m))},y_eval);
% 			lhsblk = cellfun(@(y_i){reshape(y_i(2:end,:)',[],1)},y_eval);
% 			for rhsbl_i = 1 : length(y_eval)
% 				rhsblk{rhsbl_i}(:,1:obj.n*obj.m) = cell2mat(ukrac{rhsbl_i}(1:end-1)');
% 				rhsblk{rhsbl_i}(:,obj.n*(rhsbl_i+obj.m-1)+1:obj.n*(rhsbl_i+obj.m)) = cell2mat(ac(2:size(y_eval{rhsbl_i},1))');
% 			end
% 			%--------
% 			% rhsblk = cellfun(@(y_i){zeros(numel(y_i(1:end,:)),obj.n*(length(y_eval)+1))},y_eval);
% 			% lhsblk = cellfun(@(y_i){reshape(y_i(1:end,:)',[],1)},y_eval);
% 			% for rhsbl_i = 1 : length(y_eval)
% 			%     rhsblk{rhsbl_i}(:,1:obj.n) = cell2mat(ukrac{rhsbl_i}(1:end)');
% 			%     rhsblk{rhsbl_i}(:,obj.n*rhsbl_i+1:obj.n*(rhsbl_i+1)) = cell2mat(ac(1:end)');
% 			% end
% 			%-----------
% 			rhs = cell2mat(rhsblk);
% 			lhs = cell2mat(lhsblk);
% 			% sol = rhs\lhs;
% 			sol = obj.svd_solution(lhs, rhs);
% 			% return the matrices
% 			b = reshape(sol(1:obj.n*obj.m),obj.n,obj.m);
% 			d = zeros(obj.num_obs,obj.m);
% 		end
% function x0 = initial_condition(obj, sys)
% 			% obsf  = obj.obs.obs_function;
% 			% % points = length(sys.y);
% 			% points = obj.hl_bl;
% 			% ca = obj.fullGamma(points);% Gamma
% 			% if isfield(sys, "u")
% 			% 	cab = cellfun(@(ac){ac*obj.B},ca);%
% 			% 	% Now I need a sum per sample. where each is the same size as acb:
% 			% 	% preallocate with cab
% 			% 	cabu = cab;
% 			% 	for y_k = 1 : points
% 			% 	acbuk = zeros(obj.num_obs,y_k);
% 			% 	for step = 1 : y_k
% 			% 	acbuk(:,step) = cab{y_k-step+1}*sys.u(step,:)';
% 			% 	end
% 			% 	cabu{y_k} = sum(acbuk,2);
% 			% 	end
% 			% else
% 			% 	cabu = cellfun(@(ac){zeros(obj.num_obs,1)},ca);
% 			% end
% 			%
% 			% lhs = reshape(obsf(sys.y(2:points,:))',[],1) - cell2mat(cabu(1:points-1)');
% 			% rhs = cell2mat(ca(2:points)');
% 			% x0 = obj.svd_solution(lhs, rhs);
% 			% fom this initial condition perform an optimization that finds
% 			% a better point
% 			% y = obj.C_edmd*(obj.C * x0);
% 			% x0 = pinv(obj.C)*obsf(sys.y(1,:))';
% 			x0 = pinv(obj.C)*(pinv(obj.C_edmd)*sys.y(1,:)'-obj.D*sys.u(1,:)');
% 		end
% function [cost, res, a, c] = n_cost(obj, n, U, S,Ysid,Usid)
% 			Un = U(:,1:n);
% 			Sn = S(1:n, 1:n);
% 			gam = Un*sqrtm(Sn); % Gamma
% 			gmm = gam(1:end-obj.num_obs,:); % Gamma minus
% 			gam_inv = pinv(gam); % and the inverses
% 			gmm_inv = pinv(gmm);
%       % future and past matrices
%       [y_f, u_f, w_p, yfm, ufm, wpp] = obj.fut_pst_mat(Ysid, Usid);
%
%       % zi = Yf/(Wp;Uf)
%       z_i = obj.zProj_x(y_f, [w_p; u_f]); % Uf
%       % zip = Yf-/(Wp+;Uf-)
%       if obj.det
%         zip = yfm;
%       else
%   			zip = obj.zProj_x(yfm, [wpp; ufm]);
%       end
%       % left-hand-side = (gmm_inv*zip;Yii)
%       yii = y_f(1: obj.num_obs,:);
%       % lhs = [gmm_inv*zip; ...
%       %   Ysid(obj.hl_bl*obj.num_obs+1:(obj.hl_bl+1)*obj.num_obs,:)];
%       lhs = [gmm_inv*zip; yii];
%       % right-hand-side = (gam_inv*z_i;Uf)
%       rhs = [gam_inv*z_i;u_f];
%       ac = obj.svd_solution(lhs', rhs')'; % rename U at svdDec, this is confusing
%       a = ac(1:n,1:n);
%       c = ac(n+1:n+obj.num_obs,1:n);
%       % Residuals
%       res = lhs - ac*rhs;
%       cost = sum(abs(res),'all');
%       % take the first row as the sum of state residuals
%       % res = [sum(sol(1:n,:));sol(n+1:end,:)];
%     end
% function [cost, res, a, c] = n_cost(obj, n, U, S,Ysid,Usid)
% 			Un = U(:,1:n);
% 			Sn = S(1:n, 1:n);
% 			gam = Un*sqrtm(Sn); % Gamma
% 			% gmm = gam(1:end-obj.num_obs,:); % Gamma minus
% 			gam_inv = pinv(gam); % and the inverses
% 			% gmm_inv = pinv(gmm);
%       % future and past matrices
%       [y_f, u_f, w_p, ~, ~, ~] = obj.fut_pst_mat(Ysid, Usid);
%
%       % zi = Yf/(Wp;Uf)
%       if obj.det
%         z_i = y_f;
%       else
%         z_i = obj.zProj_x(y_f, [w_p; u_f]); % Uf
%       end
%       % The state sequence is Gamma^(-1)*z_i
%       x = gam_inv*z_i;
%       % This becomes [x_{+};yii_{+}]=[A;C|K]*[x_{-};u_{f-}]
%       [x_m, x_p, yii_p, u_fm] = obj.edmd_like_division(x, y_f, u_f);
%       % zip = Yf-/(Wp+;Uf-)
%       % if obj.det
%       %   zip = yfm;
%       % else
%   		% 	zip = obj.zProj_x(yfm, [wpp; ufm]);
%       % end
%       % left-hand-side = (gmm_inv*zip;Yii)
%       % yii = y_f(1: obj.num_obs,:);
%       % lhs = [gmm_inv*zip; ...
%       %   Ysid(obj.hl_bl*obj.num_obs+1:(obj.hl_bl+1)*obj.num_obs,:)];
%       %
%       % Concatenate the matrices for the solution
%       lhs = [x_p;yii_p];
%       rhs = [x_m;u_fm];
%       % lhs = [gmm_inv*zip; yii];
%       % % right-hand-side = (gam_inv*z_i;Uf)
%       % rhs = [gam_inv*z_i;u_f];
%       ac = obj.svd_solution(lhs', rhs')'; % rename U at svdDec, this is confusing
%       a = ac(1:n,1:n);
%       c = ac(n+1:n+obj.num_obs,1:n);
%       % Residuals
%       res = lhs - ac*rhs;
%       cost = sum(abs(res),'all');
%       % take the first row as the sum of state residuals
%       % res = [sum(sol(1:n,:));sol(n+1:end,:)];
%     end
% function hlbl = hankel_blocks(obj) % (obj, yeval)
% 	% min_samples = min(cellfun(@(x) [size(x,1)], yeval));
% 	% Stablish the number of Hankel blocks
% 	% hlbl = floor(((min_samples)/((2*obj.num_obs*numel(xeval_sys))))); % Hankel blocks to use
% 	% calculate the number of blocks based on the available data
% 	% num_hl = (min_samples+1)*length(yeval);
% 	% den_hl = (2*(obj.num_obs + length(yeval)));
% 	% hlbl = floor((num_hl / den_hl)/24); % Hankel blocks to use
% 	% hlbl = floor((min_samples + 1)/(2*(obj.num_obs + 1))); % Hankel blocks to use
% 	% hlbl = floor((min_samples + 1)/(obj.num_obs + 2));
% 	% f = 1; % fatness factor
% 	% hlbl = floor(f*(min_samples + 1)/(2*(obj.num_obs + f)));
% 	% hlbl = max(ceil((min_samples + 1)/(2*(f*obj.num_obs + 1))),2);
%   % hlbl = obj.num_obs + 3;
%   hlbl = floor(2.5*obj.num_obs);
%   % hlbl = 5;
% end
% function [b, d] = BDn(obj, Ysid, Usid)
% 	% Calculates the B and D matrices of a system.
% 	if obj.m == 0
% 	b = zeros(obj.n,1);
% 	d = zeros(obj.num_obs,1);
% 	return
% 	end
% 	gam = obj.Gamma; % Gamma
% 	gmm = gam(1:end-obj.num_obs,:); % Gamma minus
% 	gam_inv = pinv(gam); % and the inverses
% 	gmm_inv = pinv(gmm);
% 	% future and past matrices
% 	[y_f, u_f, w_p, yfm, ufm, wpp] = obj.fut_pst_mat(Ysid, Usid);
% 	yii = y_f(1: obj.num_obs,:);
% 	% zi = Yf/(Wp;Uf)
% 	z_i = obj.zProj_x(y_f, [w_p; u_f]); % Uf
% 	% zip = Yf-/(Wp+;Uf-)
% 	zip = obj.zProj_x(yfm, [wpp; ufm]);
% 	% Get the Nk matrices
% 	Nk = obj.Nk(gmm, gam_inv, gmm_inv);
% 	% Get the Uk matrices
% 	Uk = arrayfun(@(bl){u_f(bl*obj.m+1:(bl+1)*obj.m,:)},0:obj.hb-1)';
% 	% Sum over all kronecker products
% 	sum_kron = obj.sum_kron(Uk,Nk);
% 	% wee need thee Pee vectorized; vec(P) = vec([gmm_inv*zip;yii]-[A;C]*gam_inv*z_i)
%   % This becomes inconsistent with the new way of calculating the A and
%   % C matrices. Waht is the fix?
% 	vec_P = reshape([gmm_inv*zip;yii]-[obj.A;obj.C]*gam_inv*z_i,[],1);
% 	% DB = sum_kron\vec_P;
% 	DB = obj.svd_solution(vec_P, sum_kron);
% 	b=reshape(DB(obj.num_obs*obj.m+1:end),[],obj.m);
% 	d=reshape(DB(1:obj.num_obs*obj.m),[],obj.m);
% end%
% function [y_f, u_f, w_p] = fut_pst_mat(obj, Ysid, Usid)
%   y_f = cellfun(@(y_i){y_i(end-obj.fb*(obj.num_obs)+1:end,:)},Ysid);
%   u_f = cellfun(@(u_i){u_i(end-obj.fb*obj.m+1:end,:)},Usid);
%   w_p = cellfun(@(ui,yi){[ui(1:obj.pb*obj.m,:);yi(1:obj.pb*obj.num_obs,:)]},Usid, Ysid);
%
% 	% yfm = cellfun(@(yi){yi((obj.nb + 1)*obj.num_obs+1:end,:)},Ysid);
% 	% ufm = cellfun(@(ui){ui((obj.nb + 1)*obj.m+1:end,:)},Usid);%Usid((obj.nb + 1)*obj.m+1:end,:);
% 	% wpp = cellfun(@(yi,ui){[ui(1:(obj.nb + 1)*obj.m,:); yi(1:(obj.nb + 1)*obj.num_obs,:)]},Ysid, Usid);%[Usid(1:(obj.nb + 1)*obj.m,:);Ysid(1:(obj.nb + 1)*obj.num_obs,:)];
% end
% function gam = fullGamma(obj, samples)
% 	gam = arrayfun(@(i){obj.C},1:samples);
% 	% iterate over ac and multiply by A fron the second blk
% 	for ac_blk = 2 : samples
% 	% This block is the prev times A
% 	gam{ac_blk} = gam{ac_blk - 1} * obj.A;
% 	end
% end
% function [a, c, k, n, ac] = ACKn(obj, Ysid, Usid)
%       % ACKn calculates the A, C, and K matrices, and the order n of the
%       % subspace.
%
%       [U, S] = obj.compute_svd(Ysid, Usid);
%
%       ns = obj.num_obs:sum(diag(S)>1e-9);
%       nv = arrayfun(@(n)obj.n_cost(n, U, S, Ysid, Usid), ns);
%       [~, midx] = min(nv);
%       % It is cheaper to test all and select the best than to optimize.
%       n = ns(midx);
%       % calculate everything that the "optimization" calculates inside. For
%       % the best of the orders n.
%       [~, res, a, c, ac] = obj.n_cost(n, U, S, Ysid, Usid);
%       SigWE = res*res'/(size(res,2)-1);
%       QQ = SigWE(1:obj.n,1:obj.n);
%       RR = SigWE(obj.n+1:obj.n+obj.num_obs,obj.n+1:obj.n+obj.num_obs);
%       SS = SigWE(1:obj.n,obj.n+1:obj.n+obj.num_obs);
%       [~,k,~,~] = idare(a',c',QQ,RR,SS); % Kalman filter ARE
%       k = k'; % Kalman gain
%     end
% function [cost, res, a, c, ac] = n_cost(obj, n, U, S,Ysid,Usid)
%      Un = U(:,1:n);
%      Sn = S(1:n, 1:n);
%      % We only need the inverse of Gamma
%      gam_inv = pinv(Un*sqrtm(Sn)); % and the inverses
%      % future and past matrices
%      [y_f, u_f, w_p] = obj.fut_pst_mat(Ysid, Usid);
%      % zi = Yf/(Wp;Uf)
%      z_i = obj.get_zi(y_f, u_f, w_p);
%      % Calculate the state sequences
%      x = cellfun(@(zi){gam_inv*zi},z_i);
%      % x = [gam_inv*z_i1,gam_inv*z_i2];
%      % The state sequence is Gamma^(-1)*z_i
%      % x = gam_inv*z_i;
%      % This becomes [x_{+};yii_{+}]=[A;C|K]*[x_{-};u_{f-}]
%      [x_p, x_f, yii_p, u_fp] = obj.edmd_like_division(x, y_f, u_f);
%      % Concatenate the matrices for the solution
%      lhs = [cell2mat(x_f);cell2mat(yii_p)];
%      rhs = [cell2mat(x_p);cell2mat(u_fp)];
%      % Get the solution
%      ac = obj.regression(lhs', rhs')'; % rename U at svdDec, this is confusing
%      a = ac(1:n,1:n);
%      c = ac(n+1:n+obj.num_obs,1:n);
%      % Residuals
%      res = lhs - ac*rhs;
%      cost = sum(abs(res),'all');
%    end
% function [b, d] = BDn(obj, Ysid, Usid)
%      % Calculates the B and D matrices of a system.
%      if obj.m == 0
%        b = zeros(obj.n,1);
%        d = zeros(obj.num_obs,1);
%        return
%      end
%      gam = obj.Gamma; % Gamma
%      gmm = gam(1:end-obj.num_obs,:); % Gamma minus
%      gam_inv = pinv(gam); % and the inverses
%      gmm_inv = pinv(gmm);
%      % future and past matrices
%      [y_f, u_f, w_p] = obj.fut_pst_mat(Ysid, Usid);
%      % get all the yii
%      yii = cellfun(@(yi){yi(1: obj.num_obs,:)},y_f);
%      % zi = Yf/(Wp;Uf)
%      z_i = obj.get_zi(y_f, u_f, w_p);
%      % zip = Yf-/(Wp+;Uf-)
%      % zip = cellfun(@(ymi,umi,wpi){obj.zProj_x(ymi,[umi; wpi])},yfm, ufm, wpp);
%      % Get the Nk matrices
%      Nk = obj.Nk(gmm, gam_inv, gmm_inv);
%      % Get the Uk matrices
%      Uk = cellfun(@(ui){arrayfun(@(bl){ui(bl*obj.m+1:(bl+1)*obj.m,:)},0:obj.fb-1)'},u_f);
%      % Sum over all kronecker products
%      sum_kron = obj.sum_kron(Uk,Nk);
%      % wee need thee Pee vectorized; vec(P) = vec([gmm_inv*zip;yii]-[A;C]*gam_inv*z_i)
%      % This becomes inconsistent with the new way of calculating the A and
%      % C matrices. Waht is the fix?
%      % Get the known state sequence
%      x = cellfun(@(zi){gam_inv*zi},z_i);
%      [x_m, x_p, yii_p, ui_p] = obj.edmd_like_division(x, yii, u_f);
%      vec_P = cellfun(@(xmi,xpi,yiii){reshape([xpi;yiii]-[obj.A;obj.C]*xmi,[],1)},x_m, x_p, yii_p);%
%      % vec_P = reshape([gmm_inv*zip;yii]-[obj.A;obj.C]*gam_inv*z_i,[],1);
%      % DB = sum_kron\vec_P;
%      DB = obj.regression(cell2mat(vec_P'), cell2mat(sum_kron'));
%      b=reshape(DB(obj.num_obs*obj.m+1:end),[],obj.m);
%      d=reshape(DB(1:obj.num_obs*obj.m),[],obj.m);
%    end
%    function Nk =  Nk(obj, gmm, gam_inv, gmm_inv)
%      % Returns an array of Nk matrices
%      % Get the La and Lc matrices for
%      La = obj.A*gam_inv;
%      Lc = obj.C*gam_inv;
%      % Upper part of Gi
%      Gup = [eye(obj.num_obs) zeros(obj.num_obs, obj.n)];
%      % Lower part of Gi
%      Gbo = [zeros((obj.fb-1)*obj.num_obs, obj.num_obs), gmm];
%      Gi = [Gup; Gbo];
%      % The Na matrices
%      % Na1 = [-La,1 M1-La,2 M2-La,3 ... M{hb-1}-La,hb)] hb=hl_bl
%      Na1 = [-La(:,1:obj.num_obs), gmm_inv - La(:,obj.num_obs+1:end)]*Gi;
%      % # All the remaining
%      Nai = arrayfun(@(bl){[gmm_inv(:,(bl-1)*obj.num_obs+1:end)-La(:,bl*obj.num_obs+1:end), zeros(obj.n,bl*obj.num_obs)]*Gi},1:obj.fb-1)';
%      % Get them into one cell
%      Na = [{Na1};Nai];
%      % for the Ncs...
%      Nc1 = [eye(obj.num_obs)-Lc(:,1:obj.num_obs), -Lc(:,obj.num_obs+1:end)]*Gi;
%      % All the remaining Ncs
%      Nci = arrayfun(@(bl){[-Lc(:,(bl-1)*obj.num_obs+1:end), zeros(obj.num_obs,(bl-1)*obj.num_obs)]*Gi},2:obj.fb)';
%      Nc = [{Nc1};Nci];
%      % Reconcatenate Nci=[Nai;Nci]
%      Nk = cellfun(@(na,nc){[na;nc]},Na,Nc);
%    end
% function skp = sum_kron(uk, nk)
%      % Performs the kronecker product and the sum of all matrices for the
%      % solution of matrices B and D.
%      kron_prod = cellfun(@(ui){cellfun(@(uii,ni){kron(uii(:,1:end-1)',ni)},ui,nk)},uk);
%      skp = c>.ellfun(@(kpi){sum(cat(3, kpi{:}), 3)},kron_prod);
%    end
