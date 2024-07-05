classdef sidDecomposition < svdDecomposition
	%SIDDECOMPOSITION subsytem identification EDMD
	properties
		K % Kalman gain of the solution
		hl_bl % number of Hankel blocks
		Gamma % observability matrix ??Return this??
		num_obs %
		C_edmd % C matrix to return to space state
	end
	methods
		function obj = sidDecomposition(observable, system)
			%SIDDECOMPOSITION
			if nargin > 0
				obj.obs = observable; % saves the observable object
				% concatenation.
				[yeval_sys, u_sys] = obj.yu_eval(system);
				% input di
				if isfield(system, 'u')
					obj.m = size(system(1).u,2);
				else
					obj.m = 0;
				end
				
				% output dim i.e., dim(y)
				obj.num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
				% nunber of hankel blocks
				obj.hl_bl = obj.hankel_blocks(yeval_sys);
				% obj.hl_bl = 10;
				% Hankelize, Notice the 2*hl_bl. Y?: The instrumental variables
				% Ysid = [Yp;Yf]
				Ysid = cell2mat(cellfun(@(x) ...
					{obj.block_hankel(x', 2*obj.hl_bl)},yeval_sys)');
				% Usid=[Up;Uf]
				if isfield(system, 'u')
					Usid = cell2mat(cellfun(@(x) ...
						{obj.block_hankel(x', 2*obj.hl_bl)},u_sys)');
				else
					Usid = zeros(2*obj.hl_bl,size(Ysid,2));
				end
				
				% Matrices A and C
				[obj.A, obj.C, obj.K, obj.n] = obj.ACKn(Ysid, Usid);
				obj.Gamma = obj.recomputeGamma(obj.A, obj.C);
				% Matrices B and D
				% [obj.B, obj.D] = obj.BD(yeval_sys, u_sys);
				[obj.B, obj.D] = obj.BDn(Ysid,Usid);
				
				obj.l = obj.obs.l;
				Cedmd = obj.matrix_C; % cannot calc and slice
				obj.C_edmd = Cedmd(:,2:end);
			end
		end
		function [cost, res, a, c] = n_cost(obj, n, U, S,Ysid,Usid)
			Un = U(:,1:n);
			Sn = S(1:n, 1:n);
			gam = Un*sqrtm(Sn); % Gamma
			gmm = gam(1:end-obj.num_obs,:); % Gamma minus
			gam_inv = pinv(gam); % and the inverses
			gmm_inv = pinv(gmm);
			% future and past matrices
			[y_f, u_f, w_p, yfm, ufm, wpp] = obj.fut_pst_mat(Ysid, Usid);
			
			% zi = Yf/(Wp;Uf)
			z_i = obj.zProj_x(y_f, [w_p; u_f]); % Uf
			% zip = Yf-/(Wp+;Uf-)
			zip = obj.zProj_x(yfm, [wpp; ufm]);
			% left-hand-side = (gmm_inv*zip;Yii)
			yii = y_f(1: obj.num_obs,:);
			% lhs = [gmm_inv*zip; ...
			%   Ysid(obj.hl_bl*obj.num_obs+1:(obj.hl_bl+1)*obj.num_obs,:)];
			lhs = [gmm_inv*zip; yii];
			% right-hand-side = (gam_inv*z_i;Uf)
			rhs = [gam_inv*z_i;u_f];
			ac = obj.svd_solution(lhs', rhs')'; % rename U at svdDec, this is confusing
			a = ac(1:n,1:n);
			c = ac(n+1:n+obj.num_obs,1:n);
			% Residuals
			res = lhs - ac*rhs;
			cost = sum(abs(res),'all');
			% take the first row as the sum of state residuals
			% res = [sum(sol(1:n,:));sol(n+1:end,:)];
		end
		function [a, c, k, n] = ACKn(obj, Ysid, Usid)
			% This will optimize over the order of the system n
			% Start with the result of the SVD
			[U, S] = obj.compute_svd(Ysid, Usid);
			s = diag(S);
			% Possible values of n up to the
			% max_n = min([obj.hl_bl,sum(s>1e-4)]);
			% max_n = max([obj.hl_bl,sum(s>5e-3)]);
			% max_n = sum(s>0);
			max_n = min(height(S),5*obj.num_obs);
			ns = obj.num_obs:max_n;
			nv = arrayfun(@(n)obj.n_cost(n, U, S, Ysid, Usid), ns);
			[~, midx] = min(nv);
			% It is cheaper to test all and select the best than tod optimize.
			n = ns(midx);
			% calculate everything that the "optimization" calculates inside. For
			% the best of the orders n.
			[~, res, a, c] = obj.n_cost(n, U, S, Ysid, Usid);
			SigWE = res*res'/(size(res,2)-1);
			QQ = SigWE(1:obj.n,1:obj.n);
			RR = SigWE(obj.n+1:obj.n+obj.num_obs,obj.n+1:obj.n+obj.num_obs);
			SS = SigWE(1:obj.n,obj.n+1:obj.n+obj.num_obs);
			[~,k,~,~] = idare(a',c',QQ,RR,SS); % Kalman filter ARE
			k = k'; % Kalman gain
		end
		
		function [b, d] = BDn(obj, Ysid, Usid)
			% Calculates the B matrix of the system. I will assume that D
			% is zero, and see if I finally get a system that does not
			% diverge.
			if obj.m == 0
				b = [];
				d = [];
				return
			end
			gam = obj.Gamma; % Gamma
			gmm = gam(1:end-obj.num_obs,:); % Gamma minus
			gam_inv = pinv(gam); % and the inverses
			gmm_inv = pinv(gmm);
			% future and past matrices
			[y_f, u_f, w_p, yfm, ufm, wpp] = obj.fut_pst_mat(Ysid, Usid);
			yii = y_f(1: obj.num_obs,:);
			% zi = Yf/(Wp;Uf)
			z_i = obj.zProj_x(y_f, [w_p; u_f]); % Uf
			% zip = Yf-/(Wp+;Uf-)
			zip = obj.zProj_x(yfm, [wpp; ufm]);
			% Get the Nk matrices
			Nk = obj.Nk(gmm, gam_inv, gmm_inv);
			% Get the Uk matrices
			Uk = arrayfun(@(bl){u_f(bl*obj.m+1:(bl+1)*obj.m,:)},0:obj.hl_bl-1)';
			% Sum over all kronecker products
			sum_kron = obj.sum_kron(Uk,Nk);
			% wee need thee Pee vectorized; vec(P) = vec([gmm_inv*zip;yii]-[A;C]*gam_inv*z_i)
			vec_P = reshape([gmm_inv*zip;yii]-[obj.A;obj.C]*gam_inv*z_i,[],1);
			% DB = sum_kron\vec_P;
			DB = obj.svd_solution(vec_P, sum_kron);
			b=reshape(DB(obj.num_obs*obj.m+1:end),[],obj.m);
			d=reshape(DB(1:obj.num_obs*obj.m),[],obj.m);
		end
		function Nk =  Nk(obj, gmm, gam_inv, gmm_inv)
			% Returns an array of Nk matrices
			% Get the La and Lc matrices for
			La = obj.A*gam_inv;
			Lc = obj.C*gam_inv;
			% Upper part of Gi
			Gup = [eye(obj.num_obs) zeros(obj.num_obs, obj.n)];
			% Lower part of Gi
			Gbo = [zeros((obj.hl_bl-1)*obj.num_obs, obj.num_obs), gmm];
			Gi = [Gup; Gbo];
			% The Na matrices
			% Na1 = [-La,1 M1-La,2 M2-La,3 ... M{hb-1}-La,hb)] hb=hl_bl
			Na1 = [-La(:,1:obj.num_obs), gmm_inv - La(:,obj.num_obs+1:end)]*Gi;
			% # All the remaining
			Nai = arrayfun(@(bl){[gmm_inv(:,(bl-1)*obj.num_obs+1:end)-La(:,bl*obj.num_obs+1:end), zeros(obj.n,bl*obj.num_obs)]*Gi},1:obj.hl_bl-1)';
			% Get them into one cell
			Na = [{Na1};Nai];
			% for the Ncs...
			Nc1 = [eye(obj.num_obs)-Lc(:,1:obj.num_obs), -Lc(:,obj.num_obs+1:end)]*Gi;
			% All the remaining Ncs
			Nci = arrayfun(@(bl){[-Lc(:,(bl-1)*obj.num_obs+1:end), zeros(obj.num_obs,(bl-1)*obj.num_obs)]*Gi},2:obj.hl_bl)';
			Nc = [{Nc1};Nci];
			% Reconcatenate Nci=[Nai;Nci]
			Nk = cellfun(@(na,nc){[na;nc]},Na,Nc);
		end
		function [u, s] = compute_svd(obj, Ysid, Usid)
			% From the hankel matrices, calculate the oblique projection
			Yf = Ysid(obj.hl_bl*obj.num_obs+1:end,:);
			Uf = Usid(obj.hl_bl*obj.m+1:end,:);
			Wp = [Usid(1:obj.hl_bl*obj.m,:);Ysid(1:obj.hl_bl*obj.num_obs,:)];
			Oi = obj.obliqueP(Yf, Uf, Wp);
			% Use the moesp way Oi*oUf the orthogonal projection of Oi onto
			% the row space of the orthogonal complement of U future
			WOW = obj.zProjOx(Oi,Uf);
			
			% Next, singular value decomp
			[u,s,~] = svd(WOW);
		end
		function gam = recomputeGamma(obj, A, C)
			gam_cell = arrayfun(@(x){C},1:obj.hl_bl)';
			for block = 2 : obj.hl_bl
				gam_cell{block} = gam_cell{block-1}*A;
			end
			gam = cell2mat(gam_cell);
		end
		function gam = fullGamma(obj, samples)
			gam = arrayfun(@(i){obj.C},1:samples);
			% iterate over ac and multiply by A fron the second blk
			for ac_blk = 2 : samples
				% This block is the prev times A
				gam{ac_blk} = gam{ac_blk - 1} * obj.A;
			end
		end
		function pred = predict(obj, x0, n_points, system)
			obsf  = obj.obs.obs_function;
			if ~isfield(system,'u')
				u = arrayfun(@(x) {zeros(x,1)}, n_points);
				obj.B = zeros(obj.n,1); % Assigns this B for the execution, it is
				% not persistent.
				obj.D = zeros(obj.num_obs,1);
			else
				u = arrayfun(@(sys){sys.u},system);
			end
			% preallocate
			% apparently, the result of an arrayfun inherits the shape of the iterable
			pred = arrayfun(@(x) struct('sv', zeros(x, obj.n), ...
				'y', zeros(x, obj.obs.l)), n_points);
			
			% Simulate
			for orb = 1 : size(x0,1)
				pred(orb).y(1,:) = obj.C_edmd*obj.C*x0{orb};
				pred(orb).sv(1,:) = x0{orb}';
				for step = 2 : n_points(orb)
					if ~mod(step, inf)%ceil(obj.hl_bl/4))
						pred(orb).sv(step,:) = obj.A*pred(orb).sv(step-1,:)' + ...
							obj.B*u{orb}(step-1,:)' + ...
							obj.K*(obsf(system(orb).y(step-1,:))' - obj.C*pred(orb).sv(step-1,:)');
					else
						pred(orb).sv(step,:) = obj.A*pred(orb).sv(step-1,:)' + ...
							obj.B*u{orb}(step-1,:)';
					end
					
					pred(orb).y(step,:) = obj.C_edmd*(obj.C*pred(orb).sv(step,:)' + obj.D*u{orb}(step,:)');
				end
			end
		end
		function x0 = initial_condition(obj, sys)
			% obsf  = obj.obs.obs_function;
			% % points = length(sys.y);
			% points = obj.hl_bl;
			% ca = obj.fullGamma(points);% Gamma
			% if isfield(sys, "u")
			% 	cab = cellfun(@(ac){ac*obj.B},ca);%
			% 	% Now I need a sum per sample. where each is the same size as acb:
			% 	% preallocate with cab
			% 	cabu = cab;
			% 	for y_k = 1 : points
			% 	acbuk = zeros(obj.num_obs,y_k);
			% 	for step = 1 : y_k
			% 	acbuk(:,step) = cab{y_k-step+1}*sys.u(step,:)';
			% 	end
			% 	cabu{y_k} = sum(acbuk,2);
			% 	end
			% else
			% 	cabu = cellfun(@(ac){zeros(obj.num_obs,1)},ca);
			% end
			%
			% lhs = reshape(obsf(sys.y(2:points,:))',[],1) - cell2mat(cabu(1:points-1)');
			% rhs = cell2mat(ca(2:points)');
			% x0 = obj.svd_solution(lhs, rhs);
			% fom this initial condition perform an optimization that finds
			% a better point
			% y = obj.C_edmd*(obj.C * x0);
			% x0 = pinv(obj.C)*obsf(sys.y(1,:))';
			x0 = pinv(obj.C)*(pinv(obj.C_edmd)*sys.y(1,:)'-obj.D*sys.u(1,:)');
		end
		function pred = pred_from_test(obj, system)
			% predict all the xts for this decomposition:
			% extract the initial conditions
			x0 = arrayfun(@(x) {obj.initial_condition(x)}, system);
			% extract the numper of points per orbit
			np = arrayfun(@(x) size(x.y,1), system);
			% make the prediction
			pred = obj.predict(x0, np, system);
		end
		function err = error(obj, system)
			% error returns the error over all the provided testing
			% samples
			% Get the prediction to calculate error
			pred = obj.pred_from_test(system);
			xts_mat = cell2mat(arrayfun(@(x) x.y, system, ...
				'UniformOutput', false));
			pred_mat = cell2mat(arrayfun(@(x) x.y, pred, ...
				'UniformOutput', false));
			err = sum(abs(xts_mat-pred_mat)./(abs(xts_mat)+eps),"all")/length(xts_mat)/obj.obs.l;
		end
		function err = abs_error(obj, xts)
			% error returns the error over all the provided testing
			% samples
			% Get the prediction to calculate error
			pred = obj.pred_from_test(xts);
			xts_mat = cell2mat(arrayfun(@(x) {x.y}, xts));
			pred_mat = cell2mat(arrayfun(@(x) {x.y}, pred));
			err = sum(abs(xts_mat-pred_mat),"all")/length(xts_mat)/obj.obs.l;
		end
		function [yev, u] = yu_eval(obj, system)
			% this is different than the others because I need the
			% evaluation with the observalbes for each sample in the
			% system, not as a complete matrix
			% ytr = arrayfun(@(x) )
			obsrv = obj.obs.obs_function;
			yev = arrayfun(@(sys) {obsrv(sys.y)}, system);
			% extract the inputs into a cell
			if ~isfield(system,'u')
				u = [];
			else
				u = {system.u}';
			end
		end
		function oi = obliqueP(obj, X, Y, Z)
			% Oi = X/_{Y}(Z) the oblique projection of X, onto
			% the row space of Z along Y. For that,
			% Oi = X/_{Y}(Z) = (X/oY)*pinv(Z/oY)*Z,
			% where the (X/oY) notation is the projection of Y onto the
			% orthogonal complement of Y
			
			% Yf/oUf
			yfPOuf = obj.zProjOx(X, Y); % Y
			% Wp/oUf, Wp=[Up;Yp]
			wpPOuf_i = pinv(obj.zProjOx(Z, Y));
			% Ready for Oi
			oi = yfPOuf*wpPOuf_i*Z;
		end
		function hlbl = hankel_blocks(obj, yeval)
			min_samples = min(cellfun(@(x) [size(x,1)], yeval));
			% Stablish the number of Hankel blocks
			% hlbl = floor(((min_samples)/((2*obj.num_obs*numel(xeval_sys))))); % Hankel blocks to use
			% calculate the number of blocks based on the available data
			num_hl = (min_samples+1)*length(yeval);
			den_hl = (2*(obj.num_obs + length(yeval)));
			hlbl = floor((num_hl / den_hl)/24); % Hankel blocks to use
			% hlbl = floor((min_samples + 1)/(2*(obj.num_obs + 1))); % Hankel blocks to use
		end
		function [y_f, u_f, w_p, yfm, ufm, wpp] = fut_pst_mat(obj, Ysid, Usid)
			% y_f: y_{future} the lower half of the y hankel matrix
			y_f = Ysid(obj.hl_bl*obj.num_obs+1:end,:);
			u_f = Usid(obj.hl_bl*obj.m+1:end,:);
			w_p = [Usid(1:obj.hl_bl*obj.m,:);Ysid(1:obj.hl_bl*obj.num_obs,:)];
			
			yfm = Ysid((obj.hl_bl + 1)*obj.num_obs+1:end,:);
			ufm = Usid((obj.hl_bl + 1)*obj.m+1:end,:);
			wpp = [Usid(1:(obj.hl_bl + 1)*obj.m,:);Ysid(1:(obj.hl_bl + 1)*obj.num_obs,:)];
		end
	end
	methods (Static)
		function skp = sum_kron(uk, nk)
			% Performs the kronecker product and the sum of all matrices for the
			% solution of matrices B and D.
			kron_prod = cellfun(@(u_i,n_i){kron(u_i',n_i)},uk,nk);
			skp = zeros(size(kron_prod{1}));
			for bl = 1 : length(kron_prod)
				skp = skp + kron_prod{bl};
			end
		end
		function zpox = zProjOx(z, x)
			% zProjOx is the orthogonal projection of z onto the row pace
			% of the orthogonal complement of x
			% xpox
			% p = size(z,1);
			% [Q,R] = qr([x;z]', "econ");
			% R = R';
			% Q = Q';
			% from the LQ decomp
			zpox = z - z*x'*pinv(x*x')*x;
			% zpox=R(end-p+1:end,end-p+1:end)*Q(:,end-p+1:end)';
		end
		function zp_x = zProj_x(z, x)
			% zProj_x is the orthogonal projection of z onto the row pace
			% of x
			% zpox = (z/_x)
			% p = size(z,1);
			% [Q,R] = qr([x;z]', "econ");
			% R = R';
			% from the LQ decomp
			% zp_x=R(end-p+1:end,1:end-p)*Q(:,1:end-p)';
			% The LQ is not numerically accurate, there is no free pizza
			zp_x = z*x'*pinv(x*x')*x;
		end
		function H = block_hankel(mat, n_blk)
			% tyr matrix of data
			% n_blk_rw is the number of block rows
			[mat_vars,j_samples] = size(mat);
			n_col = j_samples - n_blk + 1;
			H = zeros(n_blk*mat_vars, n_col);
			for blk = 1 : n_blk
				H((blk - 1)*mat_vars + 1 : blk*mat_vars, :) = ...
					mat(:,blk:blk + n_col - 1);
			end
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
