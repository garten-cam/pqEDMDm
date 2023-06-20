classdef Orbits
    %Orbits generates a set of orbits in a specified \delta{}t for any ODE
    %given to it.
    % The inputs are the structure with the necessary parameters to run the
    % solver and the solver itself
    
    
    properties (GetAccess = public, SetAccess = private)
        clean_orbits % Orbits without noise
        noise_orbits % Orbits with noise
        Conv % The equilibrium or fixed point where the orbit converges
    end
    properties(Hidden) 
        ODE = @(t,X,n,u)Orbits.DuffEqODE(t,X,struct('delta',0.5,...
                                                'beta',-1,...
                                                'alpha',1),...
                                                n,...
                                                u);
        nFunc = @Orbits.noiseFunction;
        uFunc = @Orbits.inFunction;
        Solver = 'ode45';
        % Default initial condition in a [-2 2] box
        simParam = struct('ICs',4.*rand(10,2)-2,...
                          'uP',rand(2,10),... 
                          'nP',[0 0.1],...
                          'TFin',1,...
                          'NPoints',100); 
        stop = true;
    end
    methods
        function obj = Orbits(varargin)
            %Orbits Construct an instance of this class ODE,solver,ICs,
            
            if nargin > 0 
                % P is the paramet structure for the ODE to run.
                in = inputParser;
                addOptional(in,'ODE',obj.ODE)
                addOptional(in,'nFunc',obj.nFunc)
                addOptional(in,'uFunc',obj.uFunc)
                addOptional(in,'Solver',obj.Solver,@isstr)
                addOptional(in,'simParam',obj.simParam,@isstruct)
                addOptional(in,'stop',obj.stop)
                parse(in,varargin{:})
                obj.ODE = in.Results.ODE;
                obj.nFunc = in.Results.nFunc;
                obj.uFunc = in.Results.uFunc;
                obj.Solver = in.Results.Solver;
                obj.simParam = Orbits.simParamBuilder(in.Results.simParam);
                obj.stop = in.Results.stop;
            end
            obj.clean_orbits = generatorOrbits(obj);
        end
        function clean_orb = generatorOrbits(obj)
            simSpan = 0:obj.simParam.TFin/obj.simParam.NPoints:obj.simParam.TFin;
            clean_orb = struct('Time', cell(size(obj.simParam.ICs,1),1), ...
                'SV', cell(size(obj.simParam.ICs,1),1));
            
%             [sv, time, input, fixConv] = deal(cell(size(obj.simParam.ICs,1),1));
            for inicIdx = 1:size(obj.simParam.ICs,1)
                ahODE  = obj.adhocODE(inicIdx);
                ssAdHoc = @(t,X)obj.steadyState(t,X,@(t,X)ahODE(t,X,@(t,X)obj.nFunc(t,X,[0 0])));
                if obj.stop
                    odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6,...
                    'Events',ssAdHoc);
                else
                    odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
                end
                solver = str2func(obj.Solver);
                [clean_orb(inicIdx).Time, clean_orb(inicIdx).SV] = solver(@(t,X)ahODE(t,X,@(t,X)obj.nFunc(t,X,[0 0])),...
                    simSpan,...
                    obj.simParam.ICs(inicIdx,:),...
                    odeSettings);
                
%                 clean_orb(inicIdx).U = ones(size(clean_orb(inicIdx).SV,1),size(obj.simParam.uP,1)).*obj.simParam.uP(:,inicIdx)';
%                 [~, as, ~] = obj.equilibria(inicIdx);
%                 finalSV = clean_orb(inicIdx).SV(end,:);
%                 for ias = 1:size(as,1)
%                     if all(abs(finalSV-double(as(ias,:)))<1e-1)
%                         fixConv{inicIdx} = as(ias,:);
%                     end
%                 end
            end
        end
        function [eq, as, sdd] = equilibria(obj,inIdx)
            %equilibria calculates the equilibrium/fixed points of a system
            %for a specified input index
            x = sym('x',size(obj.simParam.ICs(1,:)'),'real');
            ahODE = obj.adhocODE(inIdx);
            dSV = ahODE(0,x,@(t,X)obj.nFunc(t,X,[0 0]));
            eq = cell2sym(struct2cell(solve(dSV,x))');
            Jac = matlabFunction(jacobian(dSV),'var',{x'});
            RLam = zeros(size(eq));
            for i=1:size(eq,1)
                LamEv = eig(Jac(eq(i,:)))';
                RLam(i,:)=real(LamEv);
            end
            % Saves the AS equilibrium points
            as = eq(all(RLam<0,2),:);
            % Saves the saddle equilibrium points
            sdd = eq(not(xor(any(RLam<0,2),any(RLam>0,2))),:);
        end
    end
    methods (Hidden)
        function ahODE = adhocODE(obj,varargin)
            %adhocODE creates an ODE with an specific noise and input, this
            %is necessary because the input structure can have an array of
            %input values. This is for the generation of the orbits and the
            %calculation of fixed points and their stability
            % The user specifies a single input for the orbits, there is
            % no selection of the different inputs, therefore a varargin
            % checks for an index value for the desired input. otherwise it
            % selects the only one.
            
            in = inputParser;
            addOptional(in,'index',1);
            parse(in,varargin{:})
            inIndex = in.Results.index;
          
            ahODE = @(t,X,n)obj.ODE(t,X,n,... % General noise function 
                @(t,X)obj.uFunc(t,X,obj.simParam.uP(:,inIndex)));
            
        end
    end
    methods (Static)
        function sP = simParamBuilder(simP)
            % Builder for the simulation parameters. In case there is no
            % basic information in the simulation parameters structure.
            % at least, a matrix of initial conditions, a value for the
            % initial time and a value for the final time must be provided.
            if ~isfield(simP,'ICs')
                error('non-existent field ICs in the simulation parameters structure')
            else
                sP.ICs = simP.ICs;
            end
            if ~isfield(simP,'TFin')
                error('non-existent field TFin in the simulation parameters structure')
            else
                sP.TFin = simP.TFin;
            end
            if ~isfield(simP,'NPoints')
                error('non-existent field NPoints in the simulation parameters structure')
            else
                sP.NPoints = simP.NPoints;
            end
            if ~isfield(simP,'nP')
                sP.nP = [0 0];
            else
                sP.nP = simP.nP;
            end
            if ~isfield(simP,'uP')
                sP.uP = zeros(size(sP.ICs'));
            elseif ~(size(sP.ICs,1) == size(simP.uP,2))
                error('Number of initial conditions does not match the number of inputs')
            else
                sP.uP = simP.uP;
            end
        end
        function [value,isTerminal,direction] = steadyState(t,X,ahODE)
            value = norm(ahODE(t,X))-1e-2;
            isTerminal = 1;
            direction = 0;
        end
        function [Dx] = DuffEqODE(t,X,P,n,u)
            %DuffEqODE is the default ODE for the Koopman operator
            %algorithm. It serves an an exmaple for the implementation of
            %other ODEs.
            Dx1 = X(2);
            Dx2 = -P.delta*X(2) - P.beta*X(1) - P.alpha*X(1)^3;
            Dx = [Dx1;Dx2] + n(t,X) + u(t,X);
        end
        function n = noiseFunction(~,X,nP)
            %noiseFunction Specifies a noise or disturbance in the system
%             n = zeros(size(X,2));
            n = normrnd(nP(1),nP(2),size(X)); 
        end
        function u = inFunction(~,~,uP)
            %inFunction specifies a function for the input.
            u = uP;
        end
    end
end

