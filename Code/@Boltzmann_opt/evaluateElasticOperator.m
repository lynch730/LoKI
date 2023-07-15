function evaluateElasticOperator(bz)
    % evaluateElasticOperator is in charge of the evaluation of the
    % elements of the elastic collision operator. The elements of the
    % operator are given by the following expression:
    %
    %   g_c(u) = -N*sqrt(2*e/me)*2*u^2*sigmaC(u)
    %
    % For more information see the notes about the discretisation of the
    % boltzmann equation.
    
    % definition of intermediate variables
    kb = Constant.boltzmannInEV;      % boltzmann constant in eV
    Tg = bz.workCond.gasTemperature;  % gas temperature in K
    eg = bz.energyGrid;
    N = eg.cellNumber; % number of energy cells in the energy grid
    
    % evaluation of the elements of the operator
    bz.g_c = 2.0 * eg.node.^2.0 .* bz.elasticCrossSection;
    
    % evaluation of the boundary conditions
    bz.g_c(1) = 0;
    bz.g_c(end) = 0;
    
    % loading of the matrix
    if isempty(bz.elasticMatrix)
        bz.elasticMatrix = zeros(eg.cellNumber);
    end
    
    % Factors
    fac1 = kb*Tg/eg.step;
    fac2 = (fac1 - 0.5) / eg.step;
    fac1 = (fac1 + 0.5) / eg.step;

    % evaluate diagonal elements
    bz.elasticMatrix(1:N+1:N*N) = -(bz.g_c(1:N)*fac1 + bz.g_c(2:N+1)*fac2);

    % evaluate inferior diagonal elements
    bz.elasticMatrix(2:N+1:N*N) = bz.g_c(2:N) .* fac2;

    % evaluate superior diagonal elements
    bz.elasticMatrix(N+1:N+1:N*N) = bz.g_c(2:N) .* fac1;
    
end
