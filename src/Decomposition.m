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
classdef Decomposition
    %DECOMPOSITION parent class for all the decompositions to avoid the
    %code repetition present in the current implementation

    properties % input properties
        VVFO % input property of every decomposition
    end
    properties % calculated properties
        U % Complete decomposition matrix psi(x) = U*psi(x,u)
        A % Evolution Matrix psi(x) = A*psi(x) + B*u
        B % Input matrix              this one   ^
        C % Output matrix x = C*psi(x)
        evol_function
    end

    methods
        function obj = Decomposition(observable, xtr, ytr)
            %DECOMPOSITION Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
                obj.VVFO = observable;
                obj.U = calculate_U(obj, xtr, ytr);
                % for the A, B and C properties, I have to call a method in
                % the constructor to assign them. otherwise it always
                % executes the whole get method.
                obj.A = obj.calculate_A;
                obj.B = obj.calculate_B;
                obj.C = obj.calculate_C;
                obj.evol_function = obj.calculate_evol;
            end
        end
        function u = calculate_U(obj, xtr, ytr)
            u = g(obj, xtr)\a(obj, xtr, ytr);
        end
        function g = g(obj, xtr)
            % u'=g\a
            %    ^ this g. It also modifies the VVFO for orthogonalization
            % evaluate in the unorthogonalized observable
            x_eval = obj.xy_eval(obj.VVFO, xtr);
            % QR for orthogonalization
            [~, obj.VVFO.R] = qr(x_eval, "econ");
            % evaluate x again
            x_eval = obj.xy_eval(obj.VVFO, xtr);
            % Calculate the second g
            g = (x_eval'*x_eval)*(1/size(xtr,1));
        end
        function a = a(obj, xtr, ytr) 
            % u'=g\a < this a
            x_eval = obj.xy_eval(obj.VVFO, xtr);
            y_eval = obj.xy_eval(obj.VVFO, ytr);
            a = (x_eval'*y_eval)*(1/size(xtr,1));
        end
        function a = calculate_A(obj)
            % get the A matrix from the system x = Ax
            a = obj.U(1:end-obj.VVFO.nU, 1:end-obj.VVFO.nU)';
        end
        function b = calculate_B(obj)
            if ~obj.VVFO.nU
                b = [];
            else
                b = obj.U(end-obj.VVFO.nU+1:end,1:end-obj.VVFO.nU)';
            end
        end
        function c = calculate_C(obj)
            order_one_obs = obj.VVFO.polynomial_base(:, ...
                2:obj.VVFO.nSV + 1);
            z = sym('z',[1 obj.VVFO.nSV],'real');
            sol = solve(order_one_obs' == z');
            [AA, bb] = equationsToMatrix([struct2array(sol)]);
            c = [-bb, AA, zeros(obj.VVFO.nSV, ...
                size(obj.VVFO.polynomial_base,2) - (obj.VVFO.nSV + 1))];
        end
        function evf = calculate_evol(obj)
            % variables for the function
            x = sym('x', [1 obj.VVFO.nSV + obj.VVFO.nU]);
            % function in symbolic form
            poly_prod = obj.VVFO.polynomial_base * obj.VVFO.R_trx * obj.U;
            % function in arithmetic form
            evf = matlabFunction(poly_prod,'var',{x});
        end
        function pred = predict(obj, x0, n_points, u)
            if nargin < 4
                u = [];
            end
            % preallocate
            pred = arrayfun(@(x) struct('sv', zeros(x, obj.VVFO.nSV)), n_points);
            % assign the initial condition
            for orb = 1 : size(x0,1)
                pred(orb).sv(1,:) = x0(orb,:);
            end
            for orb = 1 : size(x0,1)
                for datapoint = 2 : n_points(orb)
                    % assign a dummy x_prev. Only necessary because there
                    % can be inputs in the system
                    if isempty(u)
                        xprev = pred(orb).sv(datapoint-1,:);
                    else
                        xprev = [pred(orb).sv(datapoint - 1, :), ...
                            u(orb).u(datapoint - 1, :)];
                    end
                    % calculate the evolution
                    xpost = obj.evol_function(xprev);
                    % If there was orhtogonalization, get back to the
                    % original space of un-orthogonalized observables
                    xpost_ogn = xpost * obj.VVFO.R;
                    pred(orb).sv(datapoint,:) = obj.C * xpost_ogn';
                end
            end
        end
        function err = error(obj, xts)
            % error returns the error over all the provided testing
            % samples
            % predict all the xts for this decomposition
            % extract the initial conditions
            x0 = cell2mat(arrayfun(@(x) x.sv(1, :), xts, ...
                'UniformOutput', false));
            % extract the numper of points per orbit
            np = arrayfun(@(x) size(x.sv,1), xts);
            % preallocate
            % !!preallocate all the things!!!
            pred = obj.predict(x0, np);
            xts_mat = cell2mat(arrayfun(@(x) x.sv, xts, ...
                'UniformOutput', false));
            pred_mat = cell2mat(arrayfun(@(x) x.sv, pred, ...
                'UniformOutput', false));
            err = sum(abs(xts_mat-pred_mat)./(abs(xts_mat)+eps),"all")/length(xts_mat)/obj.VVFO.nSV;
        end
    end
    methods (Static)
        function ev = xy_eval(obs, xytr)
            % xy_eval evaluates the data with the observable or VVFO
            Psi = obs.Psi;
            ev = [ones(size(xytr,1),1).*obs.R_trx(1,1), Psi(xytr)];
        end
    end
end