classdef hermiteObservable < pqObservable
    % HERMITEOBSERVABLE Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for l variables.
		%
		% (a, b) = (-1, 1)
    % 
    % The order of the polynomials is not greater than p, and the quasi
    % norm of the indices per polonomial function is q.
    % 
    %   O = HERMITEOBSERVABLE() 
    % 
    %   O = HERMITEOBSERVABLE(p, q, l) 
    % 
    
    methods
        function o_sym = assign_poly(obj)
            %ASSIGN_POLY Assign a polynomial basis to the specific type of
            %polynomial
						orders = obj.polynomials_order;
						% Assign Hermites
						o_sym = cell2sym( ...
							arrayfun( ...
							@(var, x_i){hermiteH(orders(var, :), x_i)}, ...
							1:obj.l, sym("x", [1 obj.l],"real"))');
        end
    end
end

