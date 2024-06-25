% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2023 Camilo Garcia-Tenorio
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this
% software and associated documentation files (the "Software"), to deal in the Software
% without restriction, including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or
% substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
classdef sidDecompositionOld < pqDecomposition
    %SIDDECOMPOSITION subsytem identification EDMD
    properties
        K % Kalman gain of the solution
        order
        hl_bl
        Gamma
        C_edmd
        
    end
    methods
        function obj = sidDecompositionOld(observable, system)
            %SIDDECOMPOSITION
            obj.obs = observable; % saves the observable object
            % concatenation.
            [xeval_sys, u_sys] = obj.xu_eval(system);
            % number of blocks

            % CCA for general purpose noise
            m = size(system(1).u,2);
            num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
            samples_per_trj = cellfun(@(x) {size(x,1)},xeval_sys);
            total_samples = min([samples_per_trj{:}]);
            hl_bl = floor(((total_samples + 1)/(2*(num_obs + 1)))/(numel(xeval_sys)/2)); % Hankel blocks to use
            % save it, for posterity
            obj.hl_bl = hl_bl;

            Usid = cellfun(@(x) ...
                {obj.block_hankel(x', 2*hl_bl)},u_sys);
            Ysid = cellfun(@(x) ...
                {obj.block_hankel(x', 2*hl_bl)},xeval_sys);

            UY = [cell2mat(Usid');cell2mat(Ysid')];
            [Q,R] = qr(UY',"econ");
            R=R'; % is the same as A = LQ';

            % Calculate L_u matrices LLL = [Lup Luf Lyp]
            % [Lup Luf Lyp] = R_[5:6,1:4]*pinv(R_[1:4,1:4])
            LLL = R(end - num_obs*hl_bl + 1:end, ... %last num_o * hl_bl i.e., Yp
                1: 2*m*hl_bl + num_obs*hl_bl ... 2*l*hlbl + num_obs*hl_bl i.e., [Up Uf Up]
                ) / ...
                (R(1:2*m*hl_bl + num_obs*hl_bl, ... 2*l*hlbl + num_obs*hl_bl i.e., [Up Yf Up]
                1:2*m*hl_bl + num_obs*hl_bl)); % 2*l*hlbl + num_obs*hl_bl i.e., [Up Yf Up]

            % First oblique projection
            % Oi = to the obliqur projection of Yf along Uf onto the row
            % space of [Up;Yp]
            % Oi = (Lup*L_[1:1,1:4] + Lyp*R_[4:4],[1:4])Q'[1:4]
            Oi = ( ...
                LLL( ...
                :,1:m*hl_bl) * ...
                R(1:m*hl_bl, 1:2*m*hl_bl + num_obs*hl_bl) + ... L_[1:1,1:4]
                LLL(:, end - num_obs*hl_bl + 1 : end) * ...
                R(2*m*hl_bl+1:2*m*hl_bl + num_obs*hl_bl, ...
                1 : 2*m*hl_bl + num_obs*hl_bl)) * ...
                Q(:, 1 : 2*m*hl_bl + num_obs*hl_bl)';
            % Second oblique projection
            % Oim meaning O_{i-1} the oblique projection of Yf- along Uf-
            % onto the row space of [Up+;Yp+]
            % First we need
            % LLLp meaning [Lup+ Luf+ Lyp+] = R_[6:1,1:5]*pinv(R_[1:5,1:5])
            LLLp = R(end - num_obs*(hl_bl-1) + 1 : end, ...
                1: 2*m*hl_bl + num_obs*(hl_bl + 1)) * ...
                pinv(R(1:2*m*hl_bl + num_obs*(hl_bl + 1), ...
                1:2*m*hl_bl + num_obs*(hl_bl + 1)));
            % Now Oim = (Lup+ * R_[1:2,1:5] + Lyp+ * R[4:5,1:5]) * Q'[1:5]
            Oim = (LLLp(:, 1:m*(hl_bl+1)) * ...
                R(1 : m*(hl_bl+1), ...
                1 : 2*m*hl_bl + num_obs*(hl_bl + 1)) + ...
                LLLp(:, end - num_obs*(hl_bl+1) + 1 : end) * ...
                R(2*m*hl_bl + 1 : 2*m*hl_bl + num_obs*(hl_bl+1), ...
                1:2*m*hl_bl + num_obs*(hl_bl + 1))) * ...
                Q(:, 1 : 2*m*hl_bl + num_obs*(hl_bl+1))';
            % From the projections, get the past and future state
            % First, get the svd that later gives Gamma
            % svd Lup * R_[1:1,1:4] + Lyp * R_[4:4,1:4]
            [U,S,~] =  svd(LLL( ...
                :,1:m*hl_bl) * ...
                R(1:m*hl_bl, 1:2*m*hl_bl + num_obs*hl_bl) + ...
                LLL(:, end - num_obs*hl_bl + 1 : end) * ...
                R(2*m*hl_bl+1:2*m*hl_bl + num_obs*hl_bl, ...
                1 : 2*m*hl_bl + num_obs*hl_bl),'econ');
            % Observability matrix Gamma
            obj.order = num_obs*hl_bl;
            % obj.order = sum(diag(S)/max(S,[],"all")>0.1);%size(S,1); % size depends on the effective rank
            U = U(:,1:obj.order);
            S = S(1:obj.order,1:obj.order);
            % V = V(:,1:obj.order);
            obj.Gamma = U*sqrtm(S);

            Xpst = (pinv(obj.Gamma)*Oi); Xpstt = Xpst';
            Xfut = pinv(obj.Gamma(1:end-num_obs,:))*Oim; Xfutt = Xfut';

            XfutY = [Xfut;UY(2*m*hl_bl+num_obs*hl_bl+1: ...
                2*m*hl_bl+num_obs*(hl_bl+1),:)];
            XpstU = [Xpst;UY(m*hl_bl+1:m*(hl_bl+1),:)];

            % use again the svd
            [Usol, Ssol, Vsol] = svd(XpstU');

            s = diag(Ssol);
            r = obj.effective_rank(s, XpstU');

            Ur = Usol(:,1:r);
            Sr = diag(s(1:r));
            Vr = Vsol(:,1:r);
            Dr = Sr\Ur'*XfutY';

            ABCD = (Vr*Dr)';



            % ABCD = futX_pstY/pstX_pstU; % System matrices
            obj.A = ABCD(1:obj.order,1:obj.order); %obj.A = A;
            obj.B = ABCD(1:obj.order,obj.order+1:obj.order+m); %obj.B = B;
            obj.C = ABCD(obj.order+1:obj.order+num_obs,1:obj.order); % the C for sid.
            obj.D = ABCD(obj.order+1:obj.order+num_obs,obj.order+1:obj.order+m); %obj.D = D;

            % Get the state error
            W = Xfut-obj.A*Xpst-obj.B*UY(m*hl_bl+1:m*(hl_bl+1),:);
            % Get the output error
            Yest = obj.C*Xpst-obj.D*UY(m*hl_bl+1:m*(hl_bl+1),:);
            % Get the corresponding output from the training set
            Y_hlbl = UY(2*m*hl_bl+num_obs*hl_bl+1: ...
                2*m*hl_bl+num_obs*(hl_bl+1),:);
            % YYYt = YYY';

            E = Y_hlbl - Yest;
            % Get the residuals
            SigWE = [W;E]*[W;E]'/(size(W,2)-1);
            QQ = SigWE(1:obj.order,1:obj.order);
            RR = SigWE(obj.order+1:obj.order+num_obs,obj.order+1:obj.order+num_obs);
            SS = SigWE(1:obj.order,obj.order+1:obj.order+num_obs);
            [P,L,G,Rept] = idare(obj.A',obj.C',QQ,RR,SS); % Kalman filter ARE
            obj.K = L'; % Kalman gain

            obj.l = obj.obs.l;
            Cedmd = obj.matrix_C;
            obj.C_edmd = Cedmd(:,2:end);
            % The orther method
            % sig = dlyap(A,QQ);
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
                    if ~mod(step, 1)%ceil(obj.hl_bl/4))
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
            num_obs = size(obj.obs.polynomials_order, 2);
            % Number of inputs
            l = size(obj.D,2);
            % The number of points to use in the approx of the ic is the same
            % as the number of Hankel blocks
            % We need the observability matrix Gamma times B
            GB = obj.Gamma*obj.B;
            hl = kron(eye(obj.hl_bl),obj.D);
            for blk_col = 0 : obj.hl_bl - 1
                hl((1 + blk_col)*num_obs + 1:end, blk_col*l + 1:l*(blk_col + 1)) = ...
                    GB(1:end-(1+blk_col)*num_obs,:);
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
