classdef pqObservable
    % pqObservable Creates an instace of a monomials based function with a
    % p-q-quasi norm reduction for l variables.
    % 
    % The order of the polynomials is not greater than p, and the quasi
    % norm of the indices per polonomial function is q.
    % 
    % Even though it is not recomendable, this observable will return a
    % monomial basis. All other typers of polynomial will inherit from this
    % class to make the corresponding implementation
    % 
    % O = pqObservable() Monomial based function for the default class
    % values
    % 
    % O = pqObservable(p, q, l)
    % 
    
    properties
        p (1,1) {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(p,10)} = 2  % Max polynomial order
        q (1,1) {mustBePositive} = 1 % q quasi-norm
        l (1,1) {mustBePositive, mustBeInteger} = 2 % Number of variables to observe    
    end

    methods
        function obj = pqObservable(p, q, l)
            % pqObservable Construct an instance of an observable. 
            if nargin > 0 
                obj.p = p;
                obj.q = q;
                obj.l = l;
                % obj.huge_p_matrix = obj.check_p_matrix(hpm, obj.n, obj.p);
            end
        end
        function p_o = polynomials_order(obj)
            pm = obj.p_matrix;
            if ~isinf(pm)
                p_o = pm( : , vecnorm(pm, obj.q) <= obj.p);
            else
                p_o = obj.elementwise_orders( ...
                    obj.l, obj.p, obj.q);
            end
            % Sort the thing lexicografically
            [~, sorting] = sort(vecnorm(p_o, obj.p));
            p_o = p_o(:, sorting);
        end
        function poly_b = poly_base(obj)
            % Array of symbolic variables: variables of the function
            xsym = sym("x", [1 obj.l],"real");
            % Get the order of the polynomials
            orders = obj.polynomials_order;
            % Preallocate array of dummy variables for the output
            o_sym = sym(ones(size(orders)));
            % Loop over the state variables to assign the polynomial
            o_sym = obj.assign_poly(o_sym, orders, xsym);
            poly_b= prod(o_sym,1);
        end
        function o_func = obs_function(obj)
            xsym = sym("x", [1 obj.l],"real");
            poly_b = obj.poly_base;
            o_func = matlabFunction(poly_b,'var',{xsym});
        end
        function o_sym = assign_poly(obj, o_sym, orders, xsym)
            %assign_poly Assign a polynomial basis to the specifiv type of
            %polynomial
            for var = 1 : obj.l
                o_sym(var,:) = (xsym(var).*o_sym(var, :)).^orders(var, :);
            end
        end
        % Check if two polynimials are eaqual based on their matrix of
        % orders.
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
 
        function hpm = p_matrix(obj)
            if  (obj.p^(obj.l)-1>9.2234e+18)
                error('Imposible to calculate observable, maximum number of iteration in elementwise calculation will be exceedeed')
            elseif (obj.l)*obj.p(1)^(obj.l) >= 1e10
                hpm = inf;
                warning('huge_pMatrix is too big for vector comparison, this is going to take while... Actually a for... you know what I mean.')
            else
                hpm = flip(dec2base(1:(obj.p(1)+1)^(obj.l) - 1, obj.p(1) + 1) - '0',2)';
            end
        end
    end
    methods (Static)
        function p_o = elementwise_orders(l, p, q)
            p_o = zeros(l,1);
            apsicount = 2;
            for i = 1 : (p+1)^(l) - 1 % This can be done in chunks
                test = flip(dec2base(i, p(1) + 1,l) - '0')';
                if vecnorm(test, q) <= p
                    p_o(:,apsicount) = test;
                    apsicount = apsicount + 1; % for indexing
                end
            end
        end
    end
end