classdef legendreObservable < pqObservable
    % legedreObservable Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for 'l' variables.
		%
		% (a, b) = (-1, 1)
    % 
    % The order of the polynomials is not greater than 'p'p, and the quasi
    % norm of the indices per polynomial function is less than 'q'.
    % 
    %   O = legendreObservable() Laguerre based function for the default class
    % values
    % 
    %   O = legendreObservable(p, q, l)
    % 
   
    methods
        function o_sym = assign_poly(obj)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
						orders = obj.polynomials_order;
						% Assign Legendre
						o_sym = cell2sym( ...
							arrayfun( ...
							@(var, x_i){legendreP(orders(var, :), x_i)}, ...
							1:obj.l, sym("x", [1 obj.l],"real"))');
        end
    end
end