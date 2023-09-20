% Author - Camilo Garcia-Tenorio Ph.D.
% MIT License

% Copyright (c) 2022 Camilo Garcia-Tenorio
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
classdef pqEDMD
    %PQEDMD calculates the pqEDMD for a particular system
    % with or without inputs
    
    properties % Inputs
        % I should not save this structure as a property, It adds
        % overhead to the class because it just keeps a copy of the
        % class arguments. The same happens with the other inputs.
        % This information should be stored in the VVFO object.
        %%%%%%% Removed System from properties
        % system struct {systemStructureCheck(system)} = struct('Time', cell(1,1), 'SV', cell(1,1))
        %%%%%%%
        p (1,:) {mustBePositive, mustBeInteger, mustBeLessThanOrEqual(p,10)} = 3
        q (1,:) {mustBePositive} = 1.2
        % [User fits whatever they want to fit] tr_ts struct {tr_ts_check(tr_ts)} = struct('tr_index', 1, 'ts_index', 2)
        polynomial (1,:) char {mustBeMember(polynomial,...
            {'Hermite',...
            'Legendre',...
            'Laguerre',...
            'ChebyshevT',...
            'ChebyshevU',...
            'Gegenbauer',...
            'Jacobi'})} = 'Legendre'
        polyParam (1,:) {mustBePositive, mustBeInteger} = [1 1]
        % [Also responsibility of the user] normalization = false
        method (1,:) char {mustBeMember(method,{'OLS',... Ordinary Least Squares
            '',... Normal Decomposition
            'maxLike' ... Maximum lilekihood decomposition
            })} = ''
    end
    properties % Calculated properties that help in the calculation
        best_pqEDMD % Best solution
    end
    properties (Dependent, Hidden) % Hidden because it is only on demand
        best_pqDecomposition
    end
    methods
        function obj = pqEDMD(system, varargin)
            if nargin > 0
                % name-value attributes for constructor
                in = inputParser;
                % addRequired(in,'system');
                addOptional(in, 'p', obj.p);
                addOptional(in, 'q', obj.q);
                addOptional(in, 'polynomial', obj.polynomial);
                addOptional(in, 'polyParam', obj.polyParam);
                addOptional(in, 'method', obj.method);
                
                parse(in, system, varargin{:});
                
                % obj.system = in.Results.system;
                obj.p = in.Results.p;
                obj.q = in.Results.q;
                % obj.tr_ts = in.Results.tr_ts;
                obj.polynomial = in.Results.polynomial;
                obj.polyParam = in.Results.polyParam;
                % obj.normalization = in.Results.normalization;
                obj.method = in.Results.method;
                
                % For the real thing, after dealing with the inputs...
                % check for polynomial parameters
                polycheck(obj.polynomial,obj.polyParam);
                % Use a fit function like the ML way. Do not calculate this
                % here. The user must create the object and then use the
                % fit function. I no longer need a "system" attribute.
                % It also has the benefit of not needing any training and
                % testing indexes. the user trains with whatever they want
                % to train with; the same for testing.
                % obj.pqEDMD_array = get_pqEDMD_array(obj);
                % obj.best_pqEDMD = get_best_pqEDMD(obj);
            end
        end
    end
    methods
        % function bestest = get_best_pqEDMD(obj)
        %     bestest = obj.pqEDMD_array(obj.best_pqDecomposition);
        % end
        function array_pqedmd = fit(obj, system)
            % Canged to fit, to be in line with the ML way of doing things
            vvfos = getVVFOs(obj, system);
            [X_train,Y_train] = pqEDMD.snapshots(system);
            % To avoid a huge switch case, I am using the str2func trick.
            % just be consistent with the naming of any decomposition class
            % or subclass.
            % Preallocate an array of decompositions
            [cell_pqedmd{1:numel(vvfos)}] = deal(str2func([obj.method,'Decomposition']));
            for pq_solution = 1 : numel(vvfos)
                cell_pqedmd{pq_solution} = cell_pqedmd{pq_solution}...
                    (vvfos(pq_solution), X_train, Y_train);
            end
            array_pqedmd = [cell_pqedmd{:}]; % turn into array
        end
        function vvfos = getVVFOs(obj, system)
            % creates an array of vvfos to evaluate the later.
            % Gets the number of state variables
            nSV = size(system(1).sv,2);
            % Gets the number of inputs.
            if isfield(system,'u')
                nU = size(system(1).u,2);
            else
                nU = 0;
            end
            % Create an emty vvfo array with the maximum value of
            % ontries, i.e., numel p * numel q. Al of these are not
            % distincst, therefore the loop ahead will deal with these
            % cases.
            unique_vvfos_index = 1;
            % Prealocates a suffieciently large vvfo array
            vvfos(numel(obj.p) * numel(obj.q)) = pqVVFO();
            % main loop over p and q values
            for p_iteration = 1 : numel(obj.p)
                % For a specific p value, all of the huge_pMatrices
                % are the same
                hpm = pqEDMD.huge_pMatrix(nSV,obj.p(p_iteration));
                for q_iteration = 1 : numel(obj.q)
                    % create a vvfo for evaluation
                    vvfo_eval = pqVVFO(obj.p(p_iteration),obj.q(q_iteration),...
                        'nSV', nSV,...
                        'nU', nU,...
                        'polynomial', obj.polynomial,...
                        'huge_pMatrix', hpm,...
                        'polyParam', obj.polyParam);
                    if unique_vvfos_index > 1
                        % Compare the evaluated vvfo with the previous
                        % ones to discard the equivalent ones.
                        if ~pqEDMD.compare_vvfo_eval(vvfo_eval, vvfos, unique_vvfos_index)
                            vvfos(unique_vvfos_index) = vvfo_eval;
                            unique_vvfos_index = unique_vvfos_index + 1;
                        end
                    else
                        vvfos(unique_vvfos_index) = vvfo_eval;
                        unique_vvfos_index = unique_vvfos_index + 1;
                    end
                end
            end
            % Delete unnecassary vvfos
            vvfos(unique_vvfos_index:end)=[];
        end
        % function best_pq = get.best_pqDecomposition(obj)
        %     [~,best_pq] = min([obj.pqEDMD_array.error]);
        % end
    end
    methods (Static)
        
        function hpm = huge_pMatrix(nSV, p_value)
            if (nSV)*p_value^(nSV) >= 1e8
                hpm=inf;
                warning('hugePmatrix is too big, reverting to elementwise computation... Go grab a coffee!')
            else
                if p_value == 1
                    hpm = eye(nSV);
                    hpm = [zeros(nSV,1) hpm];
                else
                    hpm = flip(dec2base(0:(p_value(1)+1)^(nSV) - 1, p_value(1)+1) - '0',2)';
                end
            end
        end
        function is_equal = compare_vvfo_eval(vvfo_eval, vvfos, unique_number)
            is_equal = 0;
            % If the evaluated one is equal to a previous
            % one, it is similar to the ones that came
            % just before. Therefore, the comparison
            % starts at the last unique element.
            last2first_vvfo_iterator = unique_number - 1;
            while ~is_equal && (last2first_vvfo_iterator > 0)
                if vvfos(last2first_vvfo_iterator) == vvfo_eval
                    is_equal = 1;
                else
                    last2first_vvfo_iterator = last2first_vvfo_iterator - 1;
                end
            end
        end
        function [xtr, ytr] = snapshots(system)
            % Number of trajectories from the system to populate the
            % snapshots
            orbits_number = numel(system);
            % First, store the snapshots in a cell
            [xtr_cell, ytr_cell] = deal(cell(orbits_number, 1));
            for trj = 1 : orbits_number
                % Extract the appropriate datpoints for x and y
                xtr_cell{trj} = system(trj).sv(1:end-2,:);
                ytr_cell{trj} = system(trj).sv(2:end-1,:);
                if isfield(system,'u')
                    xtr_cell{trj} = [xtr_cell{trj}, system(trj).u(1:end-2,:)];
                    ytr_cell{trj} = [ytr_cell{trj}, system(trj).u(2:end-1,:)];
                end
            end
            % Turn cells into a matrix
            xtr = cell2mat(xtr_cell);
            ytr = cell2mat(ytr_cell);
        end
    end
end

% Function to check if the system structure has the necessary fields in the structure
function polycheck(polynomial, polyParam)
if strcmp(polynomial, 'Gegenbauer') && ~isequal(size(polyParam),[1 1])
    eidType = 'pqEDMD_data:notValidPolynomialParameter';
    msgType = 'Gegenbauer Polynomials need a polyParam value.';
    throwAsCaller(MException(eidType,msgType))
end
if strcmp(polynomial, 'Jacobi') && ~isequal(size(polyParam), [1 2])
    eidTypeJ = 'pqEDMD_data:notValidPolynomialParameters';
    msgTypeJ = 'Jacobi Plonomials need a 1x2 parameter vector.';
    throwAsCaller(MException(eidTypeJ,msgTypeJ))
end
end
