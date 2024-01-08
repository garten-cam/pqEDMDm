classdef laguerreObservable < pqObservable
    % laguerreObservable Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for 'n' variables and an 'a' 
    % paramter value.
    % 
    % The order of the polynomials is not greater than 'p', and the quasi
    % norm of the indices per polynomial function is less than 'q'.
    % 
    % O = laguerreObservable() Laguerre based function for the default class
    % values
    % 
    % O = laguerreObservable(p, q, n, huge_p_matrix)
    % 
    % O = laguerreObservable( _ , huge_p_matrix) Provides the matrix of maximum
    % orders to avoid recalculation for several instances of the class with
    % the same p-value. Provide an empty array for the in-class calculation
    % e.g.,
    % O = laguerreObservable(p, q, n_vars, [])
properties
        a (1,1) {mustBeReal} = 1 % parameter a for the laguerreP polynomial
    end
    methods
        function obj = laguerreObservable(p, q, n, a)
            % laguerreObservable Construct an instance of a laguerre observable. 
            % 
            % This constructor is necessary because of the additinal 'a'
            % parameter of the method.
            if nargin == 0
                super_args = {};
            else
                % obj.a_par = a_par;
                super_args = {p, q, n};
            end
            obj@pqObservable(super_args{:});
            if nargin > 0
                obj.a = a;
            end
        end
        function o_sym = assign_poly(obj, o_sym, orders, xsym)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
            for var = 1 : obj.l
                o_sym(var,:) = laguerreL(orders(var, :), obj.a, xsym(var));
            end
        end
    end
end