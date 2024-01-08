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
classdef pqDecomposition
  %DECOMPOSITION parent class for all the decompositions.
  %
  % Returns a set of matrices characterizing a system based on the EDMD
  % algorithm.
  %
  %
  
  properties % input properties
    obs % input property of every decomposition. The observable
  end
  
  properties % calculated properties
    A % Evolution Matrix psi(x) = A*psi(x) + B*u
    B % Input matrix              this one   ^
    C % Output matrix y = C*psi(x) + D*u
    D % Input Output matrix
  end
  
  properties (Hidden)
    l % Number of variables
    m % Number of inputs
    n % Number of states in this case, number of observables
  end
  methods
    function obj = pqDecomposition(observable, system)
      %DECOMPOSITION Construct an instance of this class
      %   Detailed explanation goes here
      if nargin > 0
        
        obj.obs = observable; % saves the observable object
        
        [obs_pst, obs_fut] = obj.y_snapshots(system);
        
        if isfield(system,'u')
          obj.m = size(system(1).u, 2);
          [u_pst, u_fut] = obj.u_snapshots(system);
          obs_pst = [obs_pst u_pst];
          obs_fut = [obs_fut u_fut];
        else
          obj.m = 0;
        end
        
        % from the snapshots, get the regression
        g = obs_pst'*obs_pst;
        a = obs_pst'*obs_fut;
        U = g\a; % This is the transpose of U
        
        obj.l = obj.obs.l;
        obj.n = size(obj.obs.polynomials_order,2) + 1;
        
        % for the A, B and C properties, I have to call a method in
        % the constructor to assign them. otherwise it always
        % executes the whole get method.
        obj.A = obj.matrix_A(U);
        obj.B = obj.matrix_B(U);
        obj.C = obj.matrix_C;
        obj.D = zeros(obj.l,1);
      end
    end
    function [y_obs_pst, y_obs_fut] = y_snapshots(obj, system)
      obsrv = obj.obs.obs_function;
      y_obs_pst = cell2mat(arrayfun(@(trj) ...
        {[ones(size(trj.y,1)-2,1) obsrv(trj.y(1:end-2,:))]}, system));
      y_obs_fut = cell2mat(arrayfun(@(trj) ...
        {[ones(size(trj.y,1)-2,1) obsrv(trj.y(2:end-1,:))]}, system));
    end
    
    function a = matrix_A(obj, U)
      % get the A matrix from the system x = Ax
      a = U(1:end-obj.m, 1:end-obj.m)';
    end
    function b = matrix_B(obj, U)
      if ~obj.m
        b = zeros(obj.n,1);
      else
        b = U(end-obj.m+1:end,1:end-obj.m)';
      end
    end
    function c = matrix_C(obj)
      orders = obj.obs.poly_base;
      order_one_obs = orders(:, 1:obj.l );
      z = sym('z',[1 obj.l],'real');
      sol = solve(order_one_obs' == z');
      [AA, bb] = equationsToMatrix([struct2array(sol)],z);
      c = double([-bb, AA, zeros(obj.obs.l, ...
        size(obj.obs.poly_base,2) - obj.l)]);
    end
    function pred = predict(obj, y0, n_points, u)
      if nargin < 4
        u = arrayfun(@(x) {zeros(x,1)}, n_points);
      end
      % preallocate
      pred = arrayfun(@(x) struct('y', zeros(x, obj.obs.l)), n_points);
      % assign the initial condition
      
      obsf = obj.obs.obs_function;
      for orb = 1 : size(y0,1) % For all the initial conditions
        pred(orb).y(1,:) = y0(orb,:);
        for step = 2 : n_points(orb)
          x_prev = [1 obsf(pred(orb).y(step-1,:))];
          x_post = obj.A*x_prev' + obj.B*u{orb}(step-1,:)';
          pred(orb).y(step,:) = (obj.C*x_post + obj.D*u{orb}(step-1,:)')';
        end
      end
    end
    function err = error(obj, yts)
      % error returns the error over all the provided testing
      % samples
      % Get the prediction to calculate error
      pred = obj.pred_from_test(yts);
      % e = (y-hat(y))/y , normalized by the number of points, num of outs
      % and num of testing trajectories
      err = sum(arrayfun(@(y,hy) ... y in the training set and hat(y) in the prediction
        sum(abs(y.y-hy.y)./(abs(y.y)+eps),"all")/length(y.y)/obj.obs.l,...
        yts,pred))/numel(yts);
    end
    function err = abs_error(obj, yts)
      % error returns the error over all the provided testing
      % samples
      % Get the prediction to calculate error
      pred = obj.pred_from_test(yts);
      % e = (y-hat(y)) , normalized by the number of points, num of outs
      % and num of testing trajectories
      err = sum(arrayfun(@(y,hy) ... y in the training set and hat(y) in the prediction
        sum(abs(y.y-hy.y),"all")/length(y.y)/obj.obs.l,...
        yts,pred))/numel(yts);
    end
    function pred = pred_from_test(obj, xts)
      % predict all the xts for this decomposition:
      % extract the initial conditions
      y0 = cell2mat(arrayfun(@(x) x.y(1, :), xts, ...
        'UniformOutput', false));
      % extract the numper of points per orbit
      np = arrayfun(@(x) size(x.y,1), xts);
      % preallocate
      % I have to deal with the inputs. This was wrong....
      if isfield(xts, "u")
        pred = obj.predict(y0, np, {xts.u});
      else
        pred = obj.predict(y0, np);
      end
    end
    function spectrum(obj)
      figure()
      e = eig(obj.A);
      scatter(real(e),imag(e))
      hold on
      plot(sin(0:pi/100:2*pi),cos(0:pi/100:2*pi))
    end
  end
  methods (Static)
    function [u_pst, u_fut] = u_snapshots(system)
      u_pst = cell2mat(arrayfun(@(trj) ...
        {trj.u(1:end-2,:)}, system));
      u_fut = cell2mat(arrayfun(@(trj) ...
        {trj.u(2:end-1,:)}, system));
    end
    
  end
end
