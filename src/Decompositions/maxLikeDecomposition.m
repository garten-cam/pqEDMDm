classdef maxLikeDecomposition < Decomposition
    % Apply maximum likelihood decomposition
    properties
        Q
    end
    methods
        function obj = maxLikeDecomposition(observable, system)
            obj@Decomposition(observable, system)
            % Must recalculate Q because the parent constructor is
            % able to use it, but not store it.... ????
            [~, xtr] = obj.snapshots(system);
            obj.Q = q_matrix(obj, xtr);
        end
        function u = calculate_U(obj, system)
            % Initail condition for the optimization
            u0 = calculate_U@Decomposition(obj, system);
            [xtr, ytr] = obj.snapshots(system);
            % x and y evaluation
            x_eval = obj.xy_eval(obj.VVFO, xtr);
            % Now y
            y_eval = obj.xy_eval(obj.VVFO, ytr);
            
            % Measurement Covariance matrix
            Q_in_f = q_matrix(obj, xtr);

            options = optimset('display','iter', ...
                'Algorithm','quasi-newton');
            % Optimize
            u = zeros(size(u0));
            u(1,1) = 1;
            for u_column = 2 : size(obj.VVFO.polynomial_base,2)
                u(:,u_column) = fminunc(@(u_col)obj.cost_function(u_col, ... parameters to identify
                    Q_in_f, ... measurement covariance matrix
                    Q_in_f(u_column, u_column), ... y column covariance
                    x_eval, ... Full x_eval matrix
                    y_eval(:,u_column)), ... % The y column to approximate
                    u0(:,u_column),options);% -abs(u(:,u_column))*1000, abs(u(:,u_column))*1000,options); % used the u from OLS as starting point.
                    % x0         , lb                   , ub
            end
        end
        % function q = get.Q(obj)
        %     q = obj.Q_;
        % end
        function q = q_matrix(obj, xtr)
            x_eval = obj.xy_eval(obj.VVFO, xtr);
            x_mean = ones(1,length(x_eval))*x_eval./length(x_eval);
            meas_mean = x_eval-ones(length(x_eval),1)*x_mean;
            q = meas_mean'*meas_mean./(length(x_eval)-1);
        end
    end
    methods (Static)
        function j = cost_function(u_col, Q, sigma, x_eval, y_eval_col)
            j = sum((y_eval_col - x_eval*u_col).^2./(sigma^2 + u_col'*Q*u_col)); % x+ = Ax + Bu  ; y = Cx
        end
    end 
end