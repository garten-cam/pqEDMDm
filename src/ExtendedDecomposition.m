classdef ExtendedDecomposition
    %EXTENDEDDECOMPOSITION makes the decomposition of a 
    
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
        function obj = ExtendedDecomposition(VVFO, Xtr, Ytr, Xts, Yts, center, scale)
            %EXTENDEDDECOMPOSITION 
            % Returns the decomposition according to the observables
            % and snapshots
            if nargin > 0
                obj.VVFO = VVFO;
                obj.center = center;
                obj.scale = scale;
                obj.U = getU(obj,Xtr,Ytr);
                obj.A = getA(obj);
                obj.B = getB(obj);
                obj.h_psi = get_h_psi(obj);
                obj.C = getC(obj, Xtr);
                obj.error = get_error(obj, Xts, Yts);
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

            u = g\a; % Regularized Least squares
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
            error = sum(abs(ysim-yedmd)./(abs(ysim)+eps),"all")/length(ysim)/obj.VVFO.nSV;
        end
        function obs = get.observability(obj)
            obs = length(obj.A) - rank(obsv(obj.A, obj.C));
        end
    end
end

