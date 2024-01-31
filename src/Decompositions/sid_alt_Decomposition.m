classdef sid_alt_Decomposition < pqDecomposition
    %SIDDECOMPOSITION subsytem identification EDMD
    properties
        K % Kalman gain of the solution
        order
        hl_bl
         Gamma
    end
    properties (Hidden)
     
        num_obs
       
        C_edmd

    end

    methods
        function obj = sid_alt_Decomposition(observable, system)
            %SIDDECOMPOSITION
            obj.obs = observable; % saves the observable object
            % concatenation.
            [xeval_sys, u_sys] = obj.xu_eval(system);
            % input dim
            obj.m = size(system(1).u,2);
            % output dim i.e., dim(y)
            obj.num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
            samples_per_trj = cellfun(@(x) {size(x,1)},xeval_sys);
            total_samples = min([samples_per_trj{:}]);
            obj.hl_bl = floor(((total_samples + 1)/(2*(obj.num_obs + 1)))/(numel(xeval_sys)/2)); % Hankel blocks to use
            
            % Hankelize
            % Usid=[Up;Uf]
            Usid = cell2mat(cellfun(@(x) ...
                {obj.block_hankel(x', 2*obj.hl_bl)},u_sys)');
            % Ysid = [Yp;Yf] 
            Ysid = cell2mat(cellfun(@(x) ...
                {obj.block_hankel(x', 2*obj.hl_bl)},xeval_sys)');
            % From the hankel matrices, calculate the oblique projection
            Oi = obj.obliqueP(Usid,Ysid);
            % Use the moesp way Oi*oUf the orthogonal projection of Oi onto
            % the row space of the orthogonal complement of U future
            WOW = obj.zProjOx(Oi,Usid(obj.hl_bl*obj.m+1:end,:));

            % Next, singular value decomp
            [U,S,~] = svd(WOW);
            obj.order = sum(diag(S)/max(S,[],"all")>0.1);
            U = U(:,1:obj.order);
            S = S(1:obj.order,1:obj.order);

            % Determine the first Gamma matrices, Gamma and Gamma_m 
            obj.Gamma = U*sqrtm(S);
            
            [AC, BD, res] = obj.linear_sol_ACBD(Usid,Ysid);
            obj.A = AC(1:obj.order,1:obj.order);
            obj.C = AC(obj.order+1:obj.order+obj.num_obs,1:obj.order);
            obj = obj.recomputeGamma();
            obj.D = BD(1:obj.num_obs,:);
            obj.B = BD(obj.num_obs+1:obj.num_obs+obj.order,:);
            SigWE = res*res'/(size(res,2)-1);
            QQ = SigWE(1:obj.order,1:obj.order);
            RR = SigWE(obj.order+1:obj.order+obj.num_obs,obj.order+1:obj.order+obj.num_obs);
            SS = SigWE(1:obj.order,obj.order+1:obj.order+obj.num_obs);
            [P,L,G,Rept] = idare(obj.A',obj.C',QQ,RR,SS); % Kalman filter ARE
            obj.K = L'; % Kalman gain

            obj.l = obj.obs.l;
            Cedmd = obj.matrix_C;
            obj.C_edmd = Cedmd(:,2:end);


        end
        function obj = recomputeGamma(obj)
            obj.Gamma(1:obj.num_obs,:) = obj.C;
            for block = 2 : obj.hl_bl
                obj.Gamma((block-1)*obj.num_obs+1:block*obj.num_obs,:) = ...
                    obj.Gamma((block-2)*obj.num_obs+1:(block-1)*obj.num_obs,:)*obj.A;
            end
        end
        function [ac, bd, res] = linear_sol_ACBD(obj, Usid, Ysid)
            gmm = obj.Gamma(1:end-obj.num_obs,:);
            % And their inverses
            gam_inv = pinv(obj.Gamma);
            gmm_inv = pinv(gmm);
            % Linear equations for A C
            % z_i = Yf/(Wp;Uf)
            z_i = obj.zProj_x( ...
                Ysid(obj.hl_bl*obj.num_obs+1:end,:), ... Yf/
                [Usid(1:obj.hl_bl*obj.m,:); ... Up
                 Ysid(1:obj.hl_bl*obj.num_obs,:); ... Yp
                 Usid(obj.hl_bl*obj.m+1:end,:)]); % Uf
            % zip = Yf-/(Wp+;Uf-)
            zip = obj.zProj_x( ...
                Ysid((obj.hl_bl+1)*obj.num_obs+1:end,:), ... Yf-/
                [Usid(1:(obj.hl_bl+1)*obj.m,:); ... Up+
                 Ysid(1:(obj.hl_bl+1)*obj.num_obs,:); ... Yp+
                 Usid((obj.hl_bl+1)*obj.m+1:end,:)]); % Uf-
            % left-hand-side = (gmm_inv*zip;Yii)
            lhs = [gmm_inv*zip; ...
                Ysid(obj.hl_bl*obj.num_obs+1:(obj.hl_bl+1)*obj.num_obs,:)];
            % right-hand-side = (gam_inv*z_i;Uf)
            rhs = [gam_inv*z_i;Usid(obj.hl_bl*obj.m+1:end,:)];
            ac = lhs/rhs;
            A = ac(1:obj.order,1:obj.order);
            C = ac(obj.order+1:obj.order+obj.num_obs,1:obj.order);
            % Residuals
            res = lhs - ac*rhs;

            % Now, for B and D. just follow the alg from the book
            P = lhs - [A;C]*rhs(1:obj.order,:);
            P = P(:,1:2*obj.m*obj.hl_bl);
            Q = Usid(obj.hl_bl*obj.m+1:end,1:2*obj.m*obj.hl_bl);
            L1 = A*gam_inv;
            L2 = C*gam_inv;

            M = [zeros(obj.order,obj.num_obs),gam_inv];
            X = [eye(obj.num_obs), zeros(obj.num_obs,obj.order); ...
                 zeros(obj.num_obs*(obj.hl_bl-1),obj.num_obs), gmm];

            totm = 0;
            for k = 1 : obj.hl_bl
                N = [M(:,(k-1)*obj.num_obs+1:obj.num_obs*obj.hl_bl) - ...
                    L1(:,(k-1)*obj.num_obs+1:obj.num_obs*obj.hl_bl), ...
                    zeros(obj.order,(k-1)*obj.num_obs); ...
                    -L2(:,(k-1)*obj.num_obs+1:obj.num_obs*obj.hl_bl), ...
                    zeros(obj.num_obs,(k-1)*obj.num_obs)];
                if k==1
                    N(obj.order+1:obj.order+obj.num_obs,1:obj.num_obs) = ...
                        eye(obj.num_obs) + ...
                        N(obj.order+1:obj.order+obj.num_obs,1:obj.num_obs);
                end
                N = N*X;
                totm = totm + kron(Q((k-1)*obj.m+1:k*obj.m,:)',N);
            end
            P = P(:);
            bd_vec = totm\P;
            bd = reshape(bd_vec,obj.order+obj.num_obs,obj.m);
        end
        function oi = obliqueP(obj, Usid,Ysid)
             % Oi = Yf/_{Uf}(Wp) the oblique projection of Y future, onto
            % the row space of Wp along U future. For that,
            % Oi = Yf/_{Uf}(Wp) = (Yf/oUf)*pinv(Wp/oUf)*Wp
            
            % Yf/oUf
            yfPOuf = obj.zProjOx( ...
                Ysid(obj.hl_bl*obj.num_obs+1:end,:), ...
                Usid(obj.hl_bl*obj.m+1:end,:));
            % Wp/oUf, Wp=[Up;Yp]
            wpPOuf_i = pinv(obj.zProjOx( ...
                [Usid(1:obj.hl_bl*obj.m,:);Ysid(1:obj.hl_bl*obj.num_obs,:)], ...
                Usid(obj.hl_bl*obj.m+1:end,:)));
            % Ready for Oi 
            oi = yfPOuf*wpPOuf_i*[Usid(1:obj.hl_bl*obj.m,:);Ysid(1:obj.hl_bl*obj.num_obs,:)]; 
        end
        function pred = predict(obj, x0, n_points, system)
            obsf  = obj.obs.obs_function;
            if ~isfield(system,'u')
                u = arrayfun(@(x) {zeros(x,1)}, n_points);
            end
            % preallocate
            pred = arrayfun(@(x) struct('sv', zeros(x, obj.order), ...
                'y', zeros(x, obj.obs.l)), n_points);

            % Simulate
            for orb = 1 : size(x0,1)
                pred(orb).y(1,:) = obj.C_edmd*obj.C*x0{orb};
                pred(orb).sv(1,:) = x0{orb}';
                for step = 2 : n_points(orb)
                    if ~mod(step, ceil(obj.hl_bl/2))%ceil(obj.hl_bl/4))
                        pred(orb).sv(step,:) = obj.A*pred(orb).sv(step-1,:)' + ...
                            obj.B*system(orb).u(step-1,:)' + ...
                            obj.K*(obsf(system(orb).y(step-1,:))' - (obj.C*pred(orb).sv(step-1,:)' + obj.D*system(orb).u(step-1,:)'));
                    else
                        pred(orb).sv(step,:) = obj.A*pred(orb).sv(step-1,:)' + ...
                            obj.B*system(orb).u(step-1,:)';
                    end

                    pred(orb).y(step,:) = obj.C_edmd*(obj.C*pred(orb).sv(step,:)' + obj.D*system(orb).u(step,:)');
                end
            end
        end
        function hl = block_toeplitz(obj)
            % Number of observables
            % obj.num_obs = size(obj.obs.polynomials_order, 2);
            % Number of inputs
            l = size(obj.D,2);
            % The number of points to use in the approx of the ic is the same
            % as the number of Hankel blocks
            % We need the observability matrix Gamma times B
            GB = obj.Gamma*obj.B;
            hl = kron(eye(obj.hl_bl),obj.D);
            for blk_col = 0 : obj.hl_bl - 1
                hl((1 + blk_col)*obj.num_obs + 1:end, blk_col*l + 1:l*(blk_col + 1)) = ...
                    GB(1:end-(1+blk_col)*obj.num_obs,:);
            end
        end
        function x0 = initial_condition(obj, sys)
            obsf  = obj.obs.obs_function;
            % Calculates an estimate of the initial condition from a sequence
            % of outputs 'y' and their corresponding sequence of inputs.
            % y = Gamma * x0 + Hd * u

            % Get the Hl matrix. The toeplitz command works well for vectors,
            % not matrices
            Hl = obj.block_toeplitz;
            % reshapre the first hl_bl inputs and observed outputs
            U = reshape(sys.u(1:obj.hl_bl,:)',[],1);
            Y = reshape(obsf(sys.y(1:obj.hl_bl,:))',[],1);

            x0 = obj.Gamma\(Y - Hl*U);
        end
        function pred = pred_from_test(obj, system)
            % predict all the xts for this decomposition:
            % extract the initial conditions
            x0 = arrayfun(@(x) {obj.initial_condition(x)}, system);
            % extract the numper of points per orbit
            np = arrayfun(@(x) size(x.y,1), system);
            % preallocate
            % I have to deal with the inputs. This was wrong....
            if isfield(system, "u")
                pred = obj.predict(x0, np, system);
            else
                pred = obj.predict(x0, np);
            end
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
            xts_mat = cell2mat(arrayfun(@(x) x.sv, xts, ...
                'UniformOutput', false));
            pred_mat = cell2mat(arrayfun(@(x) x.y, pred, ...
                'UniformOutput', false));
            err = sum(abs(xts_mat-pred_mat),"all")/length(xts_mat)/obj.VVFO.nSV;
        end
        function [ytr, utr] = xu_eval(obj, system)
            % this is different than the others because I need the
            % evaluation with the observalbes for each sample in the
            % system, not as a complete matrix
            % ytr = arrayfun(@(x) )
            obsrv = obj.obs.obs_function;
            ytr = arrayfun(@(sys) {obsrv(sys.y)}, system);
            % extract the inputs into a cell
            if ~isfield(system,'u')
                utr = [];
            else
                utr = {system.u}';
            end
        end
        
    end
    methods (Static)
        
        function zpox = zProjOx(z, x)
            % zProjOx is the orthogonal projection of z onto the row pace
            % of the orthogonal complement of x
            p = size(z,1);
            [Q,R] = qr([x;z]', "econ");
            R = R';
            % from the LQ decomp
            zpox=R(end-p+1:end,end-p+1:end)*Q(:,end-p+1:end)';
        end
        function zpox = zProj_x(z, x)
            % zProjOx is the orthogonal projection of z onto the row pace
            % of the orthogonal complement of x
            p = size(z,1);
            [Q,R] = qr([x;z]', "econ");
            R = R';
            % from the LQ decomp
            zpox=R(end-p+1:end,1:end-p)*Q(:,1:end-p)';
        end
        function H = block_hankel(mat, n_blk)
            % tyr matrix of data
            % n_blk_rw is the number of block rows
            [mat_vars,j_samples] = size(mat);
            n_col = j_samples - 2*n_blk + 1;
            H = zeros(n_blk*mat_vars, n_col);
            for blk = 1 : n_blk
                H((blk - 1)*mat_vars + 1 : blk*mat_vars, :) = ...
                    mat(:,blk:blk + n_col - 1);
            end
        end
        function r = effective_rank(s, xeval_d)
            r = 1;
            while (r < size(xeval_d,2) && ...
                    s(r+1) > max(size(xeval_d))*eps*s(1))
                r = r + 1;
            end
        end
    end % ends static methods
end