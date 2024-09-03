classdef sidOlsDecomposition < sidDecomposition
  %SIDOLSDECOMPOSITION calculates the sid id with the first algorithm of
  %Overschee
  methods
    function obj = sidOlsDecomposition(fb, pb, observable, system)
      %SIDOLSDECOMPOSITION
      if nargin > 0
        obj.obs = observable; 
        [yeval_sys, u_sys] = obj.yu_eval(system);
        if isfield(system, 'u')
          obj.m = size(system(1).u,2);
          obj.unforced = false;
          obj.pb = max(1,pb); % Hankel blocks for the past
          % if forced, pb at least 1
        else
          obj.m = 0;
          obj.unforced = true;
        end
        obj.num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
        obj.fb = fb; % Hankel blocks for the future
        if ~pb
          obj.det = true;
        else
          obj.det = false;
          obj.pb = pb; % for the past
        end
        obj.l = obj.obs.l; % Number of state variables

        % Hankelize
        % Ysid = [Yp;Yf]
        Ysid = cellfun(@(sys){obj.block_hankel(sys', obj.fb, obj.pb)},yeval_sys)'; %%% the same 2hb
        % Usid=[Up;Uf]
        if obj.unforced
          % create all the empty matrices
          Usid = cellfun(@(sys){zeros(0,width(sys))},Ysid);
        else
          Usid = cellfun(@(x) ...
            {obj.block_hankel(x', obj.fb, obj.pb)},u_sys)';%%%% 2*obj.hl_bl
        end
        % Get the system
        [obj.A, obj.B, obj.C, obj.D, obj.K, obj.n] = obj.ABCDKn(Ysid, Usid);
        % back to the state
        Cedmd = obj.matrix_C; % cannot calc and slice
        obj.C_edmd = Cedmd(:,2:end);
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
      % calculate everything that the "optimization" calculates inside. For
      % the best of the orders n.
      [~, res, a, b, c, d] = obj.n_cost(n, U, S, Ysid, Usid);
      SigWE = res*res'/(size(res,2)-1);
      QQ = SigWE(1:obj.n,1:obj.n);
      RR = SigWE(obj.n+1:obj.n+obj.num_obs,obj.n+1:obj.n+obj.num_obs);
      SS = SigWE(1:obj.n,obj.n+1:obj.n+obj.num_obs);
      [~,k,~,~] = idare(a',c',QQ,RR,SS); % Kalman filter ARE
      k = k'; % Kalman gain
      if obj.m == 0
        b = zeros(n,1);
        d = zeros(obj.num_obs,1);
      end
    end

    function [cost, res, a, b, c, d] = n_cost(obj, n, U, S, Ysid, Usid)
      Un = U(:,1:n);
      Sn = S(1:n, 1:n);
      % We only need the inverse of Gamma
      gam_inv = pinv(Un*sqrtm(Sn)); % and the inverses
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
      % resab = ab_lhs - ab*ab_rhs;
      % Solve for C
      % yiif = Cxf
      c_lhs = cell2mat(yiip);
      c_rhs = cell2mat(xf);
      c = obj.regression(c_lhs', c_rhs')';
      d = zeros(height(c),width(b));
      % resc = c_lhs - cc*c_rhs;
      res = [cell2mat(xf);cell2mat(yiip)] ...
        - [a b;c d]*ab_rhs;
      cost = sum(abs(res),'all');
    end % funcion
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
    function [x_p, x_f, yii_p, u_fp, u_ff] = edmd_like_division(obj, x, y_f, u_f)
      % edmd_like_division
      % Slice and reconcatenate
      % x past
      x_p = cellfun(@(xi){xi(:,1:end-1)},x);
      % x future
      x_f = cellfun(@(xi){xi(:,2:end)},x);
      % yii future
      yii_p = cellfun(@(yi,xi){yi(1:obj.num_obs,1:width(xi)-1)},y_f,x);
      % u in the future, that is the past of yii_f
      u_fp = cellfun(@(ui,xi){ui(1:obj.m,1:width(xi)-1)},u_f,x);
      % u in the future that is the future of u_fp
      u_ff = cellfun(@(ui,xi){ui(1:obj.m,2:width(xi))},u_f,x);
    end% Function
  end % methods
end % class