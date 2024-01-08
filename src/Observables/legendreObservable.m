classdef legendreObservable < pqObservable
    % legedreObservable Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for 'l' variables.
    % 
    % The order of the polynomials is not greater than 'p'p, and the quasi
    % norm of the indices per polynomial function is less than 'q'.
    % 
    %   O = laguerreObservable() Laguerre based function for the default class
    % values
    % 
    %   O = laguerreObservable(p, q, l)
    % 
   
    methods
        function o_sym = assign_poly(obj, o_sym, orders, xsym)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
            for var = 1 : obj.l
                o_sym(var,:) = legendreP(orders(var, :), xsym(var));
            end
        end
    end
end