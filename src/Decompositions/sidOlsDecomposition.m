  classdef sidOlsDecomposition < sidDecomposition
    %SIDOLSDECOMPOSITION Summary of this class goes here
    %   Detailed explanation goes here
  
    methods
      % function obj = sidOlsDecomposition(observable, system)
      %   %SIDOLSDECOMPOSITION
      %   if nargin > 0
      %     obj.obs = observable; % saves the observable object
      %     % concatenation.
      %     [yeval_sys, u_sys] = obj.yu_eval(system);
      %     % input dim
      %     % input di
      %     if isfield(system, 'u')
      %       obj.m = size(system(1).u,2);
      %     else
      %       obj.m = 0;
      %     end
      %     % output dim i.e., dim(y)
      %     obj.num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
      %     % nunber of hankel blocks
      %     obj.hl_bl = obj.hankel_blocks(yeval_sys);
      %     % Hankelize
      %     % Ysid = [Yp;Yf]
      %     Ysid = cell2mat(cellfun(@(x) ...
      %       {obj.block_hankel(x', 2*obj.hl_bl)},yeval_sys)');
      %     % Usid=[Up;Uf]
      %     if isfield(system, 'u')
      %       Usid = cell2mat(cellfun(@(x) ...
      %         {obj.block_hankel(x', 2*obj.hl_bl)},u_sys)');
      %     else
      %       Usid = zeros(2*obj.hl_bl,size(Ysid,2));
      %     end
      %     % Get Gamma
      %     % [obj.Gamma, obj.n] = obj.computeGamma(Usid, Ysid);
      %     [obj.A, obj.C, obj.K, obj.n] = obj.ACKn(Ysid, Usid);
      %     obj.Gamma = obj.recomputeGamma(obj.A, obj.C);
      %     [obj.B, obj.D] = obj.BD(yeval_sys, u_sys);
      %     obj.l = obj.obs.l;
      %     Cedmd = obj.matrix_C; % cannot calc and slice
      %     obj.C_edmd = Cedmd(:,2:end);
      %   end
      % end
      function [cost, res, a, c] = n_cost(obj, n, U, S, Ysid, Usid)
        Un = U(:,1:n);
        Sn = S(1:n, 1:n);
        gam = Un*sqrtm(Sn); % Gamma
        gmm = gam(1:end-obj.num_obs,:); % Gamma minus
        gam_inv = pinv(gam); % and the inverses
        gmm_inv = pinv(gmm);

        [y_f, u_f, w_p, yfm, ufm, wpp] = obj.fut_pst_mat(Ysid, Usid);

        Yii = Ysid((obj.hl_bl-1)*obj.num_obs+1:(obj.hl_bl)*obj.num_obs,:);
        Uii = Usid((obj.hl_bl-1)*obj.m+1:obj.hl_bl*obj.m,:);
  
        Oi = obj.obliqueP(y_f, u_f, w_p);
        Oip = obj.obliqueP(yfm, ufm, wpp);

        WOiW = obj.zProjOx(Oi,u_f); % this removes the effect of the input
        WOipW = obj.zProjOx(Oip,ufm); % this removes the effect of the input
  
        Xi = gam_inv * WOiW;
        Xip = gmm_inv * WOipW;
 
        lhs = [Xip;Yii];
        rhs = [Xi;Uii];
        ac = obj.svd_solution(lhs', rhs')';%
        a = ac(1:n,1:n);
        c = ac(n+1:end,1:n);
        % c = sol(obj.n+1:end,1:obj.n);
        % d = sol(obj.n+1:end,obj.n+1:end);
        % c = obj.svd_solution(Yii', Xi')';
        res = lhs - ac*rhs;
        cost = sum(abs(res),'all');
      end % funcion
    end % methods
  end % class