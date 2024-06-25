% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2023 Camilo Garcia-Tenorio
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
classdef pqEDMDm
  %PQEDMD calculates the pqEDMD for a particular system
  % with or without inputs.
  % it recieves an array of parameters p, an array of parameters q and the
  % desireded decomposition.

  properties % Inputs
    p (1,:) {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(p,10)} = 3 % Array of maximum order p of the polynomials
    q (1,:) {mustBePositive} = 1.2 % Array of q quasi-norm values for the individual polynomials.
    observable = @legendreObservable % Observable object.
    dyn_dcp = @pqDecomposition % Dynamic decomposition object.
  end
  %
  % properties % Calculated properties that help in the calculation
  %   best_pqEDMD % Best solution
  % end
  %
  % properties
  %   best_pqDecomposition
  % end
  methods
    function obj = pqEDMDm(system, varargin)
      if nargin > 0
        % name-value attributes for constructor
        in = inputParser;
        % addRequired(in,'system');
        addOptional(in, 'p', obj.p);
        addOptional(in, 'q', obj.q);
        addOptional(in, 'observable', obj.observable);
        addOptional(in, 'dyn_dcp', obj.dyn_dcp);


        parse(in, system, varargin{:});

        % obj.system = in.Results.system;
        obj.p = in.Results.p;
        obj.q = in.Results.q;
        obj.observable = in.Results.observable;
        obj.dyn_dcp = in.Results.dyn_dcp;

      end
    end
  end
  methods
    function array_pqedmd = fit(obj, system)
      % Changed to fit, to be in line with the ML way of doing things
      vvfos = get_observables(obj, system);

      % run the decomposition for each of the unique obsevables.
      array_pqedmd = arrayfun(@(vfo) obj.dyn_dcp(vfo, system),vvfos);
    end
    function vvfos = get_observables(obj, system)
      % creates an array of vvfos to evaluate the later.
      % Gets the number of state variables
      l = size(system(1).y,2); % number of variables to observe
      % p q mesh
      [pm, qm] = meshgrid(obj.p, obj.q);
      vvfos = reshape(...
        arrayfun(@(ip,iq)obj.observable(ip,iq,l),pm,qm),1,[]);
      % populate the parameters and compare with the previous ones

      uniques = zeros(numel(vvfos), 1);
      for obs = 1 : numel(vvfos)
        uniques(obs) = obj.compare(vvfos, obs);
      end
      % Delete unnecassary vvfos
      vvfos(~uniques)=[];
    end
  end
  methods (Static)
    function unique = compare(vvfos, iter)
      % compare returns a true value if the observable of the iteration
      % is not equal to the previous ones.
      unique = 1;
      last2first_vvfo_iterator = iter - 1;
      while unique && (last2first_vvfo_iterator > 0)
        if vvfos(last2first_vvfo_iterator) == vvfos(iter)
          unique = 0;
        else
          last2first_vvfo_iterator = last2first_vvfo_iterator - 1;
        end
      end
    end

  end
end