classdef sidAltDecomposition < sidDecomposition
    %SIDALTDECOMPOSITION The subspace decomp with the least squares alg
    methods
        function obj = sidAltDecomposition(observable, system)
            %SIDALTDECOMPOSITION Construct an instance of this class
            %   Detailed explanation goes here
            obj.obs = observable; % saves the observable object
            % concatenation.
            [yeval_sys, u_sys] = obj.yu_eval(system);
            % input dim
            obj.m = size(system(1).u,2);
            % output dim i.e., dim(y)
            obj.num_obs = size(obj.obs.polynomials_order, 2); % number of outputs, number of observables
            % nunber of hankel blocks
            obj.hl_bl = obj.hankel_blcks(yeval_sys);
            % Hankelize
            % Usid=[Up;Uf]
            Usid = cell2mat(cellfun(@(x) ...
                {obj.block_hankel(x', 2*obj.hl_bl)},u_sys)');
            % Ysid = [Yp;Yf]
            Ysid = cell2mat(cellfun(@(x) ...
                {obj.block_hankel(x', 2*obj.hl_bl)},yeval_sys)');
            % Get Gamma
        end
        
        function [gam, order] = computeGamma(obj, Usid, Ysid)
            Oi = obj.obliqueP(Usid, Ysid);
            % What do I do? I need to fix the python thing also...
        end
    end
end

