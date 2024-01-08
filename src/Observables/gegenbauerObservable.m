classdef gegenbauerObservable < pqObservable
    %gegenbauerObservables Creates an instace of a Gegenbauer polynomials based
    % function with a p-q-quasi norm reduction for 'l' variables and an 'a' 
    % paramter value.
    %
    % The order of the polynomials is not greater than p, and the quasi
    % norm of the indices per polynomial function is less than q.
    %
    %   O = gegenbauerObservable()
    %
    %   O = gegenbauerObservable(p, q, l, a) 
    
    properties
        a (1,1) {mustBeReal} = 1 % parameter a for the gegenbauerC polynomial
    end
    methods
        function obj = gegenbauerObservable(p, q, l, a)
            % gegenbauerObservable Construct an instance of a gegenbauer
            % based observable. This constructor is necessary because of
            % the additinal a_parameter of the method.
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
    end
    methods
        function o_sym = assign_poly(obj, o_sym, orders, xsym)
            %assign_poly Assign a polynomial basis to the specific type of
            %polynomial
            for var = 1 : obj.l
                o_sym(var,:) = gegenbauerC(orders(var, :), obj.a, xsym(var));
            end
        end
    end
end

