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
classdef rrrDecomposition < Decomposition
    %RRRDECOMPOSITION Implements the R^3 method from Zlatko Drmac
    methods
        function obj = rrrDecomposition(observable, xtr, ytr)
            obj@Decomposition(observable, xtr, ytr)
        end
        function u = calculate_U(obj, xtr, ytr)
            xeval_d = obj.xy_eval(obj.VVFO, xtr);
            yeval_d = obj.xy_eval(obj.VVFO, ytr);
            [U,S,V] = svd(xeval_d);
            s = diag(S);
            % Efffective rank r of xeval_d
            r = 1;
            while (r < size(xeval_d,2) && s(r+1) > max(size(xeval_d))*eps*s(1))
                r = r + 1;
            end
            d = U'*yeval_d;
            u = V*([d(1:r,:)./s(1:r); zeros(size(xeval_d,2)-r,r)]);
            % TODO [ ] continue with the QR steps and still get an U mat
        end
    end
end

