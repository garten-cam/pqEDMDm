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
classdef svdDecomposition < pqDecomposition
    %SVDDECOMPOSITION Implements the svd decomposition for the calculation
    %of the state transition matrix.
    methods
        function obj = svdDecomposition(observable, system)
            if nargin > 0

                obj.obs = observable; % saves the observable object

                [obs_pst, obs_fut] = obj.y_snapshots(system);

                if isfield(system,'u')
                    obj.m = size(system(1).u, 2);
                    [u_pst, u_fut] = obj.u_snapshots(system);
                    obs_pst = [obs_pst u_pst];
                    obs_fut = [obs_fut u_fut];
                else
                    obj.m = 0;
                end
                U = obj.U(obs_pst, obs_fut);
                obj.l = obj.obs.l;
                obj.n = size(obj.obs.polynomials_order,2) + 1;

                % for the A, B and C properties, I have to call a method in
                % the constructor to assign them. otherwise it always
                % executes the whole get method.
                obj.A = obj.matrix_A(U);
                obj.B = obj.matrix_B(U);
                obj.C = obj.matrix_C;
                obj.D = zeros(obj.l,obj.m);
            end
        end
        function u = U(obj, obs_pst, obs_fut)
            [Ud,S,V] = svd(obs_pst);
            s = diag(S);
            % Efffective rank r of xeval_d
            r = obj.effective_rank(s, obs_pst);
            Ur = Ud(:,1:r);
            Sr = diag(s(1:r));
            Vr = V(:,1:r);
            Dr = Sr\Ur'*obs_fut;

            % fill dr with zeros in case r<d
            % Dd = [Dr ; zeros(size(xeval_d,2)-r,size(xeval_d,2))];
            u = Vr*Dr;
        end
    end
    methods (Static)
        function r = effective_rank(s, xeval_d)
            r = 1;
            while (r < size(xeval_d,2) && ...
                    s(r+1) > max(size(xeval_d))*eps*s(1))
                r = r + 1;
            end
        end
    end
end

