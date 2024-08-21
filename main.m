function [X, Residue, Par_abs, Par_scaled, fval] = constrainOM(B, z, D, Val_real, Val_initial, Par_true, Iterpar, measure_noisy, Wi)
    % This function performs constrained optimization for a transformed system of equations.
    % Inputs:
    %   B, z, D - Matrices/vectors representing the system of equations
    %   Val_real - Exact parameter values
    %   Val_initial - Initial guess for the optimization
    %   Par_true - True parameters to be compared with
    %   Iterpar - Iteration parameters for the optimization
    %   measure_noisy - Noisy measurement values
    %   Wi - Weighting factors for the measurements (optional)
    %
    % Outputs:
    %   X - Optimized parameters
    %   Residue - Residuals of the equations
    %   Par_abs - Absolute parameter values
    %   Par_scaled - Scaled parameter values
    %   fval - Final value of the objective function

    global MAN;

    if nargin < 9
        Wi = ones(size(B, 1), 1);
    end

    % Step 1: Obtain measured displacements in either B or D
    Measured = unique([GetMeasured(B); GetMeasured(D)]);
    Nm = length(Measured);

    % Step 2: Introduce the Scaling factor
    % Initialize and assign initial values
    Scale_factor = zeros(size(z));
    assignInitialValues(Val_initial);

    % Locate the measured values in the Val_real vector
    ind_Measur = locateMeasurements(Measured, Val_real);
    
    % Assign noisy measurements to the corresponding variables
    assignNoisyMeasurements(Measured, ind_Measur, measure_noisy);

    % Calculate the scaling factors
    Scale_factor = calculateScaleFactor(z, Scale_factor);

    % Scale matrix B and evaluate D
    d_val = eval(D);
    B_unscaled = eval(B);
    B_scaled = scaleMatrix(B_unscaled, Scale_factor);

    % Step 3: Solve the optimization problem
    [Xval, Xsym, fval] = MinimConstrainOM(B_scaled, z, d_val, Iterpar, Wi);

    % Step 4: Post-process the results
    [X, Residue] = postProcessResults(Xval, Xsym, z, Val_initial, B_scaled, d_val);

    % Step 5: Compare and scale parameters
    [Par_abs, Par_scaled] = compareAndScaleParameters(X, Par_true, Xsym, Val_initial);

    % Final step: Handle shear area error adjustment
    X = adjustForShearArea(X, Xsym, MAN);

end

%% Helper Functions

function mset = GetMeasured(M)
    % This function extracts the measured displacements from matrix M
    mset = symvar(M);
    mset = mset(:);
    ind = cellfun(@(x) ismember(x(1), {'u', 'v', 'w'}), mset);
    mset(ind) = [];
end

function assignInitialValues(Val_initial)
    % Assigns initial values to the corresponding variables
    for i = 1:size(Val_initial, 1)
        eval([char(Val_initial{i, 1}) ' = ' num2str(Val_initial{i, 2}) ';']);
    end
end

function ind_Measur = locateMeasurements(Measured, Val_real)
    % Finds the location of measurements in the Val_real vector
    ind_Measur = cellfun(@(x) ~isempty(strfind(char(Measured), char(x))), Val_real(:, 1), 'UniformOutput', false);
    ind_Measur = find(cell2mat(ind_Measur));
end

function assignNoisyMeasurements(Measured, ind_Measur, measure_noisy)
    % Assigns noisy measurement values to the corresponding variables
    for ii = 1:length(Measured)
        eval([char(Measured(ii)) ' = measure_noisy(' num2str(ind_Measur(ii)) ');']);
    end
end

function Scale_factor = calculateScaleFactor(z, Scale_factor)
    % Calculates the scaling factor for each variable in z
    for i = 1:length(z)
        aux = char(z(i));
        if length(aux) == 5
            Scale_factor(i) = eval(sym(aux));
        elseif length(aux) == 10
            Scale_factor(i) = eval(sym(aux(1:5))) * eval(sym(aux(6:10)));		
        end
    end
end

function B_scaled = scaleMatrix(B_unscaled, Scale_factor)
    % Scales the matrix B with the given scaling factors
    B_scaled = ones(size(B_unscaled));
    for i = 1:size(B_unscaled, 2)
        B_scaled(:, i) = B_unscaled(:, i) * Scale_factor(i);
    end
end

function [X, Residue] = postProcessResults(Xval, Xsym, z, Val_initial, B_scaled, d_val)
    % Processes the optimization results and calculates residuals
    ind = locate(Xsym, z);
    x_val = Xval(ind);
    Residue = B_scaled * x_val - d_val;
    
    % Find single variables and sort them
    ind_single = cellfun(@(x) length(x) == 5, Xsym);
    [Xsym, order] = sort(Xsym(ind_single));
    Xval = Xval(ind_single);
    Xval = Xval(order);
    
    ind2 = locate(Val_initial(:, 1), Xsym);
    Xval = Xval .* cell2mat(Val_initial(ind2, 2));
    
    X = [Xsym, num2cell(Xval)];
end

function [Par_abs, Par_scaled] = compareAndScaleParameters(X, Par_true, Xsym, Val_initial)
    % Compares and scales the parameters
    Par_abs = Par_true;
    
    ind1 = locate(X(:, 1), Par_true(:, 1));
    Par_abs(:, 2) = X(ind1, 2);
    
    Par_scaled = Par_abs;
    Par_scaled(:, 2) = num2cell(cell2mat(Par_abs(:, 2)) ./ cell2mat(Par_true(:, 2)));
end

function X = adjustForShearArea(X, Xsym, MAN)
    % Adjusts the shear area error in the optimization results
    shearerror = 1e-16;
    v = 0.25;
    
    for i = 1:length(Xsym)
        if startsWith(Xsym{i}, 'Q')
            auxnumQ = str2double(Xsym{i}(2:end));
            Xsym{i} = strrep(Xsym{i}, 'Q', 'As');
            valQ = cell2mat(X(i, 2));
            
            for j = 1:length(Xsym)
                if startsWith(Xsym{j}, 'I')
                    auxnumI = str2double(Xsym{j}(2:end));
                    if auxnumQ == auxnumI
                        valI = cell2mat(X(j, 2));
                        auxabs = abs(valQ);
                        if auxabs < shearerror
                            X{i, 2} = 0;
                        else
                            auxl = MAN(auxnumQ);
                            L = 1;
                            auxl = eval(auxl);
                            X{i, 2} = -(12 * valI * (2 * v + 2) * (valQ - 1)) / (valQ * auxl^2);
                        end
                    end
                end
            end
        end
    end
end

function ind1 = locate(var1, var2)
    % Locates the index of var2 in var1
    ind1 = cellfun(@(x) find(ismember(var1, x)), var2, 'UniformOutput', false);
    ind_empty = cellfun(@(x) ~isempty(x), ind1);
    ind1 = cell2mat(ind1(ind_empty));
end
