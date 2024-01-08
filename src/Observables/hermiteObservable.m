classdef hermiteObservable < pqObservable
    % HermiteObservable Creates an instace of a Hermite polynomials based
    % function with a p-q-quasi norm reduction for l variables.
    % 
    % The order of the polynomials is not greater than p, and the quasi
    % norm of the indices per polonomial function is q.
    % 
    %   O = hermiteObservable() 
    % 
    %   O = hermiteObservable(p, q, l) 
    % 
    
    methods
        function o_sym = assign_poly(obj, o_sym, orders, xsym)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
            for var = 1 : obj.l
                o_sym(var,:) = hermiteH(orders(var, :), xsym(var));
            end
        end
    end
end

