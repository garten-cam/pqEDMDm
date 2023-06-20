classdef NCA_Decomposition
    %UNTITLED Neighborhood component analysis for regression
    properties % Input properties
        VVFO
        center
        scale
    end
    properties % Calculated properties
        U % Complete matrix [psi+, u+]' = U[psi, u]'
        % Separated model
        % psi+ = A*psi + B*u
        % x = h(psi) as the inverse of monomial functions
        regdata % Measurement covariance matrix. In this case is the state covariance matrix
%         A
%         B
%         h_psi
%         % for recovering the state in a linear way
%         % x = C*psi
%         C
%         % Error of the test set
        error
    end
    properties (Dependent, Hidden)
        observability
    end
    methods
        function obj = NCA_Decomposition(VVFO, Xtr, Ytr, Xts, Yts, center, scale)
            %MaxLikeDecomp performs a maximum likelihood decompostion
            %of the EDMD
            if nargin > 0
                obj.VVFO = VVFO;
                obj.center = center;
                obj.scale = scale;
                %                 [obj.regdata, obj.U] = getRegU(obj,Xtr,Ytr);
                [obj.regdata, obj.U] = getRegU(obj,Xtr,Ytr,Xts,Yts);
                %                 obj.A = getA(obj);
                %                 obj.B = getB(obj);
                %                 obj.h_psi = get_h_psi(obj);
                %                 obj.C = getC(obj, Xtr);
                %                 obj.error = get_error(obj, Xts, Yts);
            end
        end
        function [regInfo, u] = getRegU(obj, X, Y, Xt, Yt)
            %GETU produces the U matrix of the decomposition

            % Evaluates the snapshots with the Psi function of the
            % observables

            Psi = obj.VVFO.Psi;

            % Evaluate the training and testing sets
            x_eval = Psi(X);
            y_eval = Psi(Y);

            xt_eval = Psi(Xt{1});
            yt_eval = Psi(Yt{1});

            num_observables = size(obj.VVFO.polynomial_base,2) - 1;
            % Preallocate the u matrix
            u = zeros(num_observables,num_observables);
            % the first column of u is always ones
            %u(1,1) = 1;
            % preallocate fitinfo structure
            %fitinfo(num_observables) = struct();
            % get the regularized regression for each variable
            for u_column = num_observables : -1 : 1
%                 lambdavals = linspace(0,3,50)*std(y_eval(:,u_column))/length(y_eval(:,u_column));
%                 lossvals = zeros(length(lambdavals),1);
%                 for lambda = 1:length(lambdavals)
%                     nca = fsrnca(x_eval,y_eval(:,u_column),...
%                         'FitMethod','exact', ...
%                         'Solver','sgd','Verbose',0,'Lambda',lambdavals(lambda), ...
%                         'LossFunction','epsiloninsensitive');
%                     % evalueate the loss on the testing set.
%                     lossvals(lambda) = loss(nca,xt_eval,yt_eval(:,u_column),'LossFunction','mse');
%                 end
%                 [~,id_best_lambda] = min(lossvals);
                regInfo(u_column).fitdata = fsrnca(x_eval,...
                    y_eval(:,u_column),...
                    'LossFunction','epsiloninsensitive',...
                    'FitMethod','exact',...
                    'Solver','lbfgs');
                if find(regInfo(u_column).fitdata.W~=1)
                    fprintf("w IS NOT 1")
                end
                u(:,u_column) = regInfo(u_column).fitdata.FeatureWeights;
            end
        end
%         function u = getU(obj, X, Y)
%             %GETU produces the U matrix of the decomposition
%             % Evaluates the snapshots with the Psi function of the
%             % observables
%             Psi = obj.VVFO.Psi;
%             x_eval = [ones(size(X,1),1), Psi(X)];
%             y_eval = [ones(size(Y,1),1), Psi(Y)];
%             % Calculates the G and A matrices
%             g = (x_eval'*x_eval)*(1/size(X,1));
%             % This is not the A matrix of the system, I use this
%             % notation to be consuistent with the literature
%             while rcond(g) <= eps
%                 [~,~,E] = qr(g,0);
%                 % I need to guarantee that the order one polynomials
%                 % stay on the basis. Only fo this one, If I do it with
%                 % C, there is no problem.
%                 % Get the indices of the order one poly
%                 while any(E(end)==find(sum(obj.VVFO.polynomials_order)==1))
%                     E = circshift(E,1);
%                 end
%                 obj.VVFO.polynomials_order = obj.VVFO.polynomials_order(:,sort(E(1:end-1)));
%                 Psi = obj.VVFO.Psi;
%                 if all(~logical(obj.VVFO.polynomials_order(:,1)))
%                     x_eval = [ones(size(X,1),1), Psi(X)];
%                     y_eval = [ones(size(Y,1),1), Psi(Y)];
%                 else
%                     x_eval = Psi(X);
%                     y_eval = Psi(Y);
%                 end
%                 g = (x_eval'*x_eval)*(1/size(X,1));
%             end
%             a = (x_eval'*y_eval)*(1/size(Y,1));
%             % THis returns the transpose of K because it avoids
%             % numerical errors if it is considered a Koopman operator
%             % and the espectral decomposition is necessary
%             u = g\a;
%         end
%         function a = getA(obj)
%             a = obj.U(1:end-obj.VVFO.nU, 1:end-obj.VVFO.nU)';
%         end
%         function b = getB(obj)
%             if ~obj.VVFO.nU
%                 b = [];
%             else
%                 b = obj.U(end-obj.VVFO.nU+1:end,1:end-obj.VVFO.nU);
%             end
%         end
%         function h = get_h_psi(obj)
%             obs_capture_B = eye(size(obj.U));
%             obs_capture_B = obs_capture_B(1:end-obj.VVFO.nU, sum(obj.VVFO.polynomials_order)==1);
%             obs_capture_B = obs_capture_B(:,1:end-obj.VVFO.nU);
%             % according to B, extract the polynomials that capture the state
%             sym_state_observables = obj.VVFO.polynomial_base(1:end-obj.VVFO.nU)*obs_capture_B;
%             % create auxiliary variable for the inversse
%             y = sym('y',[1 length(obs_capture_B)],'Real');
%             sym_inverse = subs(symvar(sym_state_observables),...
%                 solve(sym_state_observables==y(sum(obj.VVFO.polynomials_order(:,1:end-obj.VVFO.nU))==1)));
%             %             sym_inverse =
%             h = matlabFunction(sym_inverse,'vars',{y}) ;
%         end
%         function c = getC(obj, X)
%             Psi = obj.VVFO.Psi;
%             if all(~logical(obj.VVFO.polynomials_order(:,1)))
%                 x_eval = [ones(size(X,1),1), Psi(X)];
%             else
%                 x_eval = Psi(X);
%             end
%             c = X(:,1:end-obj.VVFO.nU)'*pinv(x_eval');
%         end
%         function error = get_error(obj, X, Y)
%             % Calculate the evolution of every initial condition in
%             % the test set according to the decomposition
%             % Create the function of observables, this case does not
%             % need vector evaluation, so the constant first term is
%             % irrelevant
%             poly_base = obj.VVFO.polynomial_base;
%             Psi = matlabFunction(poly_base,'var',{symvar(poly_base)});
%             Y_edmd = cell(length(Y),1); % preallocate
%             lagMult = cell(length(Y),1); % preallocate the lagrange multiplier
%             % Evaluate initial condition
%             for orbit = 1 : length(Y)
%                 % evaluate the initial condition of each orbit
%                 Y_edmd{orbit}(1,:) = obj.h_psi(Psi(X{orbit}(1,:))*obj.U);
%                 for datapoint = 2 : length(Y{orbit})
%                     if ~obj.VVFO.nU % If there is no input, there is no need to concatenate anything for the evaulation
%                         Y_edmd{orbit}(datapoint,:) = obj.h_psi(Psi(Y_edmd{orbit}(datapoint-1,:))*obj.U);
%                     else
%                         % add the input to the state variable
%                         state_plus_input = [Y_edmd{orbit}(datapoint-1,:), Y{orbit}(datapoint-1,end-obj.VVFO.nU+1:end)];
%                         Y_edmd{orbit}(datapoint,:) = obj.h_psi(Psi(state_plus_input)*obj.U);
%                     end
%                 end
%             end
%             % unpack both cells
%             ysim = cell2mat(Y);
%             % slice ysim in case there are inputs
%             ysim = ysim(:,1:end-obj.VVFO.nU);
%             yedmd = cell2mat(Y_edmd);
%             % Get the error
%             %             error = sum(abs(ysim-yedmd)./(abs(ysim)+eps),"all")/length(ysim)/obj.VVFO.nSV;
%             error = sum(abs(ysim-yedmd),"all")/length(ysim)/obj.VVFO.nSV;
%         end
%         function obs = get.observability(obj)
%             obs = length(obj.A) - rank(obsv(obj.A, obj.C));
%         end
    end
end

% %             % Preallocate matrux U
% %             num_observables = size(obj.VVFO.polynomial_base,2);
% % %             if num_observables > 19
% % %                 stop = 1;
% % %             end
% %             u = zeros(num_observables);
% %             u(1,1) = 1;
% %             u_col_cov = zeros(num_observables,num_observables,num_observables);
% %             for u_col = 2 : size(obj.VVFO.polynomial_base,2)
% %                 [u_col_cov(:,:,u_col),~,u(:,u_col)] = hacMod(x_eval,y_eval(:,u_col),'intercept',false,'weights','PZ','Bandwidth',1,'SmallT',false);
% %             end