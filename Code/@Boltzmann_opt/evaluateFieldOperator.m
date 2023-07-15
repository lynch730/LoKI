function evaluateFieldOperator(bz, ~, ~)
    % evaluateFieldOperator is in charge of the evaluation of the elements of the field operator. The elements of
    % the operator is given by the following expression: (please note that the factor EoN^2 is not included because
    % of compatibility reasons with the simulation of field pulses, this factor is later added when doing calculations
    % with g_E and/or fieldMatrix)
    %
    %      g_E(u) = -N*sqrt(2*e/me)*EoN^2/(sigmaT(u)*(1+WoN^2/(2*e*u*sigmaT^2(u)/me)))
    %
    % For more information see the notes about the discretisation of the boltzmann equation.
    
    %% Intermediate variables
    gamma = Constant.gamma;                % gamma parameter (sqrt(2e/me))
    WoN = bz.workCond.reducedExcFreqSI;    % reduced angular exitation frequency (SI units)
    eg = bz.energyGrid;
    N = eg.cellNumber;          % number of energy cells in the energy grid

    % evaluation of the elements of the operator
    bz.g_E = eg.node .* bz.totalCrossSection .^ 2.0;
    bz.g_E = 1.0 + (WoN./gamma).^2.0 ./ bz.g_E;
    bz.g_E = eg.node ./ ( 3.0 .* bz.totalCrossSection .* bz.g_E);
    
    % evaluation of the boundary conditions
    bz.g_E(1) = 0;
    bz.g_E(end) = 0;
    
    % loading of the matrix
    if isempty(bz.fieldMatrix)
        bz.fieldMatrix = zeros(eg.cellNumber);
    end

    % Inverse Step Size Squared
    de2_inv = 1.0 ./ (eg.step .* eg.step);

    % evaluate diagonal elements
    bz.fieldMatrix(1:N+1:N*N) = - (bz.g_E(1:N)+bz.g_E(2:N+1)) .* de2_inv;

    % evaluate inferior diagonal elements
    bz.fieldMatrix(2:N+1:N*N) = bz.g_E(2:N) .* de2_inv;

    % evaluate superior diagonal elements
    bz.fieldMatrix(N+1:N+1:N*N) = bz.g_E(2:N) .* de2_inv;

end