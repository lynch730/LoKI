function evaluateCAROperator(bz)
    % evaluateCAROperator is in charge of the evaluation of the elements of
    % the continuous approximation for rotation (CAR) operator with
    % Chapman-Cowling correction. The operator is developed for gases well
    % described by the rotational cross sections proposed by Gerjuoy and
    % Stein (complete development of the operator can be found in "M A
    % Ridenti, L L Alves, V Guerra and J Amorim, Plasma Sources Sci.
    % Technol. 24 (2015) 035002 (16pp)).
    %
    % The elements of the operator are given by the following expression:
    %
    %     g_car(u) = -N*sqrt(2*e/me)*4*u*sigma0*B
    %
    % For more information see the notes about the discretisation of the
    % boltzmann equation.
    
    % definition of intermediate variables
    kb = Constant.boltzmannInEV;                  % boltzmann constant in eV
    Tg = bz.workCond.gasTemperature;       % gas temperature in K
    eg = bz.energyGrid;
    N = eg.cellNumber; % number of energy cells in the energy grid
    
    %% evaluation of the elements of the operator

    % Loop CAR gases
    sigma0B = 0;
    for gasName = bz.CARgases

        % Return gas
        gasID = Gas.find(gasName, bz.gasArray);
        gas = bz.gasArray(gasID);

        % Calculate each
        sig_tmp = gas.electricQuadrupoleMoment^2.0;
        sig_tmp = sig_tmp * gas.fraction * gas.rotationalConstant;

        % Add to sum
        sigma0B = sigma0B + sig_tmp;

    end

    % Scale
    fac0 = 15.0 * Constant.electronCharge^2.0 * Constant.bohrRadius^2.0;
    sigma0B = 8.0 * pi * sigma0B / fac0;
    bz.g_car = 4.0 * eg.node .* sigma0B;
    
    % evaluation of the boundary conditions
    bz.g_car(1) = 0;
    bz.g_car(end) = 0;
    
    % loading of the matrix
    if isempty(bz.CARMatrix)
        bz.CARMatrix = zeros(N);
    end
    
    % Diag Factors
    fac1 = kb * Tg / eg.step;
    fac2 = (fac1 - 0.5) ./ eg.step;
    fac1 = (fac1 + 0.5) ./ eg.step;
    
    % evaluate diagonal elements
    bz.CARMatrix(1:N+1:N*N) = - (bz.g_car(1:N)*fac1+bz.g_car(2:N+1)*fac2);

    % evaluate inferior diagonal elements
    bz.CARMatrix(2:N+1:N*N) = bz.g_car(2:N).*fac2;

    % evaluate superior diagonal elements
    bz.CARMatrix(N+1:N+1:N*N) = bz.g_car(2:N).*fac1;

end