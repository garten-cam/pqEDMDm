classdef laguerreObservable < pqObservable
    % laguerreObservable Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for 'l' variables and an 'a' 
    % paramter value.
		%
		% (a, b) = (0, inf)
    % 
    % The order of the polynomials is not greater than 'p', and the quasi
    % norm of the indices per polynomial function is less than 'q'.
    % 
    % O = laguerreObservable() Laguerre based function for the default class
    % values
    % 
    % O = laguerreObservable(p, q, l, a)
properties
        a (1,1) {mustBeReal} = 1 % parameter a for the laguerreP polynomial
    end
    methods
			function obj = laguerreObservable(p, q, l, a)
            % laguerreObservable Construct an instance of a laguerre observable. 
            % 
            % This constructor is necessary because of the additinal 'a'
            % parameter of the method.
            if nargin == 0
                super_args = {};
            else
                % obj.a_par = a_par;
                super_args = {p, q, l};
            end
            obj@pqObservable(super_args{:});
            if nargin > 0
                obj.a = a;
            end
        end
        function o_sym = assign_poly(obj)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
						orders = obj.polynomials_order;
						% Preallocate array of dummy variables for the output
            o_sym = cell2sym( ...
							arrayfun( ...
							@(var, x_i){laguerreL(orders(var, :), obj.a, x_i)}, ...
							1:obj.l, sym("x", [1 obj.l],"real"))');
        end
    end
end