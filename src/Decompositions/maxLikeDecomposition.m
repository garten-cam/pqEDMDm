classdef maxLikeDecomposition < svdDecomposition
    % Apply maximum likelihood decomposition
    properties
        Q
    end
    methods
				function u = regression(obj, obs_pst, obs_fut)
            % Initail condition for the optimization
            u0 = obj.regression@svdDecomposition(obs_pst, obs_fut);
            
            % Measurement Covariance matrix
            obj.Q = obs_pst;

            options = optimset('display','iter', ...
                'Algorithm','quasi-newton');
            % Optimize
            u = zeros(size(u0));
            u(1,1) = 1;
            for u_column = 2 : size(obj.obs.poly_base,2)
                u(:,u_column) = fminunc(@(u_col)obj.cost_function(u_col, ... parameters to identify
                    obj.Q, ... measurement covariance matrix
                    obj.Q(u_column, u_column), ... y column covariance
                    obs_pst, ... Full x_eval matrix
                    obs_fut(:,u_column)), ... % The y column to approximate
                    u0(:,u_column),options);% -abs(u(:,u_column))*1000, abs(u(:,u_column))*1000,options); % used the u from OLS as starting point.
                    % x0         , lb                   , ub
            end
        end
        % function q = get.Q(obj)
        %     q = obj.Q_;
        % end
        function obj = set.Q(obj, obs_pst)
            %q_matrix gets the 
            x_mean = ones(1,length(obs_pst))*obs_pst./length(obs_pst);
            meas_mean = obs_pst-ones(length(obs_pst),1)*x_mean;
            obj.Q = meas_mean'*meas_mean./(length(obs_pst)-1);
        end
    end
    methods (Static)
        function j = cost_function(u_col, Q, sigma, x_eval, y_eval_col)
            j = sum((y_eval_col - x_eval*u_col).^2./(sigma^2 + u_col'*Q*u_col)); % x+ = Ax + Bu  ; y = Cx
        end
    end 
end