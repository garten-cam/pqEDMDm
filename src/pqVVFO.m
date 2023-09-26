% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2022 Camilo Garcia-Tenorio
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
classdef pqVVFO < handle
    %PQVVFO Creates an instace of a vector valued function of
    %observables for a single p  and q values and a type of polynomial

    properties
        p (1,1) {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(p,10)}  = 1 % Max poly order
        q (1,1) {mustBePositive}  = 1 % q norm
        nSV (1,1) {mustBePositive, mustBeInteger} = 3 % Number of state variables
        nU (1,1) {mustBeNonnegative, mustBeInteger} = 0 % Number of inputs
        polynomial (1,:) char {mustBeMember(polynomial,...
            {'Hermite',...
            'Legendre',...
            'Laguerre',...
            'ChebyshevT',...
            'ChebyshevU',...
            'Gegenbauer',...
            'Jacobi'})} = 'Legendre' % Type of polynomial
        polyParam (1,:) {mustBePositive, mustBeInteger} = [1 1] % Parameter of Gengen and Jac polynomials
        huge_pMatrix (:,:) = [] % Equivalent p values have hthe same base matrix. Input to avoid recalculation
        
    end
    properties (Dependent) %
        polynomials_order % Final matrix of polynomial orders after pq reduction
        polynomial_base % Symbolic array of observable functions
        Psi % matlabFunction for the evaluation of the observables
        R  % R matrix from QR decomposition
        R_trx % R transformation matrix. just the inverse or R
    end
    properties (Hidden, Access=private)
        % This defreats the whole purpose of having a handle class and
        % dependent methods. Now I have to update the polynomial base and
        % the Psi every time there is a change in the orders. e.g.,
        % eliminating observables with a qr decomposition.
        polynomials_order_ % To avoid recalculation when calling
        polynomial_base_ % Same ^
        Psi_ % Same             ^
        R_
        R_trx_
    end

    methods
        function obj = pqVVFO(p, q, varargin)
            %PQVVFO Construct an instance of this class
            % not a matlab.System class, because I need to have an
            % array of these objects at some point.
            % This will have its own input parser method
            if nargin > 0 % Only builds the whole object if it gets
                obj.p = p;
                obj.q = q;

                in = inputParser;
                addOptional(in,'nSV',obj.nSV,@isnumeric)
                addOptional(in,'nU',obj.nU)
                addOptional(in,'polynomial',obj.polynomial)
                addOptional(in,'polyParam', obj.polyParam);
                addOptional(in,'huge_pMatrix', obj.huge_pMatrix);
                addOptional(in,'R',[])

                parse(in,varargin{:})
                obj.nSV = in.Results.nSV;
                obj.nU = in.Results.nU;
                obj.polynomial = in.Results.polynomial;
                obj.polyParam = in.Results.polyParam;
                % Check if the huge p matrix is given, else, build it.
                obj.huge_pMatrix = obj.check_pMatrix(in.Results.huge_pMatrix, obj.nSV, obj.p);
                %
                % obj.polynomials_order = [];
                % obj.polynomial_base = [];
                obj.R = in.Results.R;
                
            end
        end
        function poly_orders = get.polynomials_order(obj)
            % I deeply dislike this trick. 
            if ~isempty(obj.polynomials_order_)
                poly_orders = obj.polynomials_order_;
            else % Assign the value to the dummy before returning
                obj.polynomials_order = [];
                poly_orders = obj.polynomials_order_;
            end
        end
        function set.polynomials_order(obj, poly_order)
            if ~isempty(poly_order)
                obj.polynomials_order_ = poly_order;
            else
                if ~isinf(obj.huge_pMatrix)
                    poly_orders = obj.huge_pMatrix(:, ...
                        vecnorm(obj.huge_pMatrix,obj.q)<=obj.p);
                else
                    poly_orders = pqVVFO.elementwisePsiOrders(obj.nSV, obj.p, obj.q);
                end
                % Add the ordering
                [~,  sorting] = sort(vecnorm(poly_orders,obj.p));
                poly_orders = poly_orders(:,sorting);
                % Add the inputs
                poly_orders = [poly_orders, zeros(size(poly_orders, 1), obj.nU);...
                    zeros(obj.nU,size(poly_orders, 2)), eye(obj.nU)];
                obj.polynomials_order_ = poly_orders;
            end
        end
        function base = get.polynomial_base(obj)
            % Returns a function to evaluate the vvfo according to the
            % order indices and the number of state variables and
            % inputs.
            if ~isempty(obj.polynomial_base_)
                base = obj.polynomial_base_;
            else
                obj.polynomial_base = [];
                base = obj.polynomial_base_;
            end
        end
        function set.polynomial_base(obj, val)
            if ~isempty(val)
                obj.polynomial_base_ = val;
            else
                % First, get the base polynomial function
                base_polyF = pqVVFO.polyString2Func(obj.polynomial, obj.polyParam);
                % Create an array of symbolic variables
                xsym = sym('x',[1 (obj.nSV + obj.nU)],'real');
                % Preallocate the matrix of symbolic variables
                xPsi = sym(ones(size(obj.polynomials_order)));
                % Loop over the state variables to assign the polynomial
                % according to the order and variable
                for state_variable = 1 : obj.nSV
                    xPsi(state_variable,1:end-obj.nU) = base_polyF(obj.polynomials_order(state_variable,1:end-obj.nU), xsym(state_variable)) ;
                end
                % Add the variables for the inputs
                if obj.nU
                    xPsi(end-obj.nU+1:end, end-obj.nU+1:end) = diag(xsym(end-obj.nU+1:end))+~eye(obj.nU);
                end
                obj.polynomial_base_ = prod(xPsi,1);
            end
        end
        function psi = get.Psi(obj)
            % Temoporal copy of the symbolic basis to avoid calling
            % the get method several times
            if ~isempty(obj.Psi_)
                psi = obj.Psi_;
            else
                obj.Psi = [];
                psi = obj.Psi_;
            end
        end
        function set.Psi(obj, val)
            if ~isempty(val)
                obj.Psi_ = val;
            else
                poly_base = obj.polynomial_base*obj.R_trx;
                if all(~logical(obj.polynomials_order(:,1)))
                    % If the first column of polybase is all zeros, then there
                    % is a constant term. if that is the case, the psi
                    % matlabFunction does not allow for vectorial evaluation
                    obj.Psi_ = matlabFunction(poly_base(2:end),'var',{sym('x',[1 length(poly_base)-1])});
                else
                    obj.Psi_ = matlabFunction(poly_base,'var',{sym('x',[1 length(poly_base)])});
                end
            end
        end
        function r = get.R(obj)
            r = obj.R_;
        end
        function rtrx = get.R_trx(obj)
            rtrx = obj.R_trx_;
        end
        function set.R(obj, r_val)
            [~, n_obs] = size(obj.polynomials_order);
            if ~isempty(r_val) &&  all(size(r_val))==all([n_obs, n_obs])
                obj.R_ = r_val;
                obj.R_trx_ = inv(r_val);
            else
                [obj.R_, obj.R_trx_] = deal(eye(n_obs));
            end
        end
        function E = eq(obj1,obj2)
            if size(obj1.polynomials_order)==size(obj2.polynomials_order)
                if all(obj1.polynomials_order==obj2.polynomials_order)
                    E = 1;
                else
                    E = 0;
                end
            else
                E = 0;
            end
        end
    end
    methods(Static)
        function hpm = check_pMatrix(hpm_in_parser, nSV, p_value)
            if ~isempty(hpm_in_parser)
                hpm = hpm_in_parser;
            elseif  (p_value^(nSV)-1>9.2234e+18)
                error('Imposible to calculate observable, maximum number of iteration in elementwise calculation will be exceedeed')
            elseif (nSV)*p_value(1)^(nSV) >= 1e10
                hpm=inf;
                warning('huge_pMatrix is too big for vector comparison, reverting to element-wise calculation. Go grab a coffee!')
            else
                hpm = flip(dec2base(0:(p_value(1)+1)^(nSV) - 1, p_value(1) + 1) - '0',2)';
            end
        end
        function aPsi = elementwisePsiOrders(nSV, p_value, q_value)
            aPsi = zeros(nSV,1);
            apsicount = 2;
            for i = 1 : (p_value+1)^(nSV) - 1 % This can be done in groups. It is too much time
                test = flip(dec2base(i, p_value(1) + 1,nSV) - '0')';
                if vecnorm(test, q_value)<=p_value
                    aPsi(:,apsicount) = test;
                    apsicount = apsicount + 1; % for indexing
                end
            end
        end
        function base_polyF = polyString2Func(poly,polyParam)
            if strcmpi(poly,'Hermite')
                poly_str = 'hermiteH';
            elseif strcmpi(poly,'Legendre')
                poly_str = 'legendreP';
            elseif strcmpi(poly,'Laguerre')
                poly_str = ['@(n,x)laguerreL(n,',...
                    num2str(polyParam(1)),',x)'];
            elseif strcmpi(poly,'ChebyshevT')
                poly_str = 'chebyshevT';
            elseif strcmpi(poly,'ChebyshevU')
                poly_str = 'chebyshevU';
            elseif strcmpi(poly,'Gegenbauer')
                poly_str = ['@(n,x)gegenbauerC(n,',...
                    num2str(polyParam(1)),',x)'];
            elseif strcmpi(poly,'Jacobi')
                poly_str = ['@(n,x)jacobiP(n,',...
                    num2str(polyParam(1)),',',...
                    num2str(polyParam(2)),',x)'];
            end
            base_polyF = str2func(poly_str);
        end
    end
end

