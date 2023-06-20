classdef RegularizedDecomposition
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

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
        A
        B
        h_psi
        % for recovering the state in a linear way
        % x = C*psi
        C
        % Error of the test set
        error
    end
    properties (Dependent, Hidden)
        observability
    end
    methods
        function obj = RegularizedDecomposition(VVFO, Xtr, Ytr, Xts, Yts, center, scale)
            %MaxLikeDecomp performs a maximum likelihood decompostion
            %of the EDMD
            if nargin > 0
                obj.VVFO = VVFO;
                obj.center = center;
                obj.scale = scale;
                [obj.regdata, obj.U] = getRegU(obj,Xtr,Ytr);
                obj.A = getA(obj);
                obj.B = getB(obj);
                obj.h_psi = get_h_psi(obj);
                obj.C = getC(obj, Xtr);
                obj.error = get_error(obj, Xts, Yts);
            end
        end
        function [fitinfo, u] = getRegU(obj, X, Y)
            %GETU produces the U matrix of the decomposition

            % Evaluates the snapshots with the Psi function of the
            % observables
            
            Psi = obj.VVFO.Psi;
            if all(~logical(obj.VVFO.polynomials_order(:,1)))
                x_eval = [ones(size(X,1),1), Psi(X)];
                y_eval = [ones(size(Y,1),1), Psi(Y)];
            else
                x_eval = Psi(X);
                y_eval = Psi(Y);
            end
            num_observables = size(obj.VVFO.polynomial_base,2);
            % Preallocate the u matrix
            u = zeros(num_observables,num_observables);
            % the first column of u is always ones
            %u(1,1) = 1;
            % preallocate fitinfo structure
            %fitinfo(num_observables) = struct();
            % get the regularized regression for each variable
            for u_column = num_observables : -1 : 1
%                 [b, fitdata] = lasso(x_eval(:,2:end),y_eval(:,u_column),...
%                     "CV",10,... do it with cross validation
%                     "Intercept",true,... no intercept, we already have that
%                     "LambdaRatio",1e-5,... try very small lambdas
%                     "Standardize",false);
                [b, fitdata] = lassoglm(x_eval(:,2:end),y_eval(:,u_column),...
                    "normal",...
                    "CV",10,... do it with cross validation
                    "LambdaRatio",1e-5,... try very small lambdas
                    "DFmax",2,... Set a max of 2 
                    "Standardize",false);                
                fitinfo(u_column).fitdata = fitdata;
                % we also need to store the b matrix to retrieve the other
                % values in the feature selection process.
                fitinfo(u_column).b = b;
                indextostore = fitdata.IndexMinDeviance;
                u(1,u_column) = fitdata.Intercept(indextostore);
                u(2:end,u_column) = b(:,indextostore);
            end
        end
        function u = getU(obj, X, Y)
            %GETU produces the U matrix of the decomposition
            % Evaluates the snapshots with the Psi function of the
            % observables
            Psi = obj.VVFO.Psi;
            x_eval = [ones(size(X,1),1), Psi(X)];
            y_eval = [ones(size(Y,1),1), Psi(Y)];
            % Calculates the G and A matrices
            g = (x_eval'*x_eval)*(1/size(X,1));
            % This is not the A matrix of the system, I use this
            % notation to be consuistent with the literature
            while rcond(g) <= eps
                [~,~,E] = qr(g,0);
                % I need to guarantee that the order one polynomials
                % stay on the basis. Only fo this one, If I do it with
                % C, there is no problem.
                % Get the indices of the order one poly
                while any(E(end)==find(sum(obj.VVFO.polynomials_order)==1))
                    E = circshift(E,1);
                end
                obj.VVFO.polynomials_order = obj.VVFO.polynomials_order(:,sort(E(1:end-1)));
                Psi = obj.VVFO.Psi;
                if all(~logical(obj.VVFO.polynomials_order(:,1)))
                    x_eval = [ones(size(X,1),1), Psi(X)];
                    y_eval = [ones(size(Y,1),1), Psi(Y)];
                else
                    x_eval = Psi(X);
                    y_eval = Psi(Y);
                end
                g = (x_eval'*x_eval)*(1/size(X,1));
            end
            a = (x_eval'*y_eval)*(1/size(Y,1));
            % THis returns the transpose of K because it avoids
            % numerical errors if it is considered a Koopman operator
            % and the espectral decomposition is necessary
            u = g\a;
        end
        function a = getA(obj)
            a = obj.U(1:end-obj.VVFO.nU, 1:end-obj.VVFO.nU)';
        end
        function b = getB(obj)
            if ~obj.VVFO.nU
                b = [];
            else
                b = obj.U(end-obj.VVFO.nU+1:end,1:end-obj.VVFO.nU);
            end
        end
        function h = get_h_psi(obj)
            obs_capture_B = eye(size(obj.U));
            obs_capture_B = obs_capture_B(1:end-obj.VVFO.nU, sum(obj.VVFO.polynomials_order)==1);
            obs_capture_B = obs_capture_B(:,1:end-obj.VVFO.nU);
            % according to B, extract the polynomials that capture the state
            sym_state_observables = obj.VVFO.polynomial_base(1:end-obj.VVFO.nU)*obs_capture_B;
            % create auxiliary variable for the inversse
            y = sym('y',[1 length(obs_capture_B)],'Real');
            sym_inverse = subs({sym("x",size(sym_state_observables))},...
                solve(sym_state_observables==y(sum(obj.VVFO.polynomials_order(:,1:end-obj.VVFO.nU))==1)));
            %             sym_inverse =
            h = matlabFunction(sym_inverse,'vars',{y}) ;
        end
        function c = getC(obj, X)
            Psi = obj.VVFO.Psi;
            if all(~logical(obj.VVFO.polynomials_order(:,1)))
                x_eval = [ones(size(X,1),1), Psi(X)];
            else
                x_eval = Psi(X);
            end
            c = X(:,1:end-obj.VVFO.nU)'*pinv(x_eval');
        end
        function error = get_error(obj, X, Y)
            % Calculate the evolution of every initial condition in
            % the test set according to the decomposition
            % Create the function of observables, this case does not
            % need vector evaluation, so the constant first term is
            % irrelevant
              Y_edmd = cell(length(Y),1); % preallocate
           % This was old code that did not get updated.......
           evol_function = matlabFunction(obj.h_psi(obj.VVFO.polynomial_base*obj.U),'var',{sym("x",size(obj.VVFO.polynomial_base))});
            % Evaluate initial condition
            for orbit = 1 : length(Y)
                % evaluate the initial condition of each orbit
                Y_edmd{orbit}(1,:) = evol_function(X{orbit}(1,:));
                for datapoint = 2 : length(Y{orbit})
                    if ~obj.VVFO.nU % If there is no input, there is no need to concatenate anything for the evaulation
                        Y_edmd{orbit}(datapoint,:) = evol_function(Y_edmd{orbit}(datapoint-1,:));
                    else
                        % add the input to the state variable
                        state_plus_input = [Y_edmd{orbit}(datapoint-1,:), Y{orbit}(datapoint-1,end-obj.VVFO.nU+1:end)];
                        Y_edmd{orbit}(datapoint,:) = evol_function(state_plus_input);
                    end
                end
            end
            % unpack both cells
            ysim = cell2mat(Y);
            % slice ysim in case there are inputs
            ysim = ysim(:,1:end-obj.VVFO.nU);
            yedmd = cell2mat(Y_edmd);
            % Get the error
            %             error = sum(abs(ysim-yedmd)./(abs(ysim)+eps),"all")/length(ysim)/obj.VVFO.nSV;
            error = sum(abs(ysim-yedmd),"all")/length(ysim)/obj.VVFO.nSV;
        end
        function obs = get.observability(obj)
            obs = length(obj.A) - rank(obsv(obj.A, obj.C));
        end
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