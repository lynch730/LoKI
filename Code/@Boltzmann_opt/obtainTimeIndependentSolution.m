function obtainTimeIndependentSolution(bz)   
    % obtainTimeIndependentSolution solves the time-independent electron
    % boltzmann equation (to be used for steady-state simulations)
    
    % Local copyu
    eg = bz.energyGrid;
    isnonlinear = bz.includeNonConservativeIonization || ...
                  bz.includeNonConservativeAttachment || ...
                  bz.includeEECollisions;

    % save appropiate method for the selected non-linear algorithm
    switch bz.nonLinearAlgorithm
        case 'mixingDirectSolutions'
            nonLinearSolver = str2func('mixingDirectSolutions');
        case 'temporalIntegration'
            nonLinearSolver = str2func('temporalIntegration');
    end
    
    % check for the presence of non-linear operators in the Boltzmann
    % equation
    if isnonlinear
        auxEedf = nonLinearSolver(bz); 
    else
        auxEedf = bz.linearSolver;
    end
    
    % when the smart grid is activated the minimum number of decades of
    % decay in the eedf is ensured
    if eg.isSmart

        % Compute EEDF Fall in Decades
        decades = log10(auxEedf(1)) - log10(auxEedf(end));

        % Iterate Solution until convergence
        while decades < eg.minEedfDecay

            % increase maximum value of the energy grid
            eg.updateMaxValue(eg.node(end)*(1+eg.updateFactor));
            
            % check for the presence of non-linear operators in the
            % Boltzmann equation
            if isnonlinear
                auxEedf = nonLinearSolver(bz);
            else
                auxEedf = bz.linearSolver;
            end
            
            % check decades of decay
            decades = log10(auxEedf(1)) - log10(auxEedf(end));

        end

        % Iterate Solution until convergence
        while decades > eg.maxEedfDecay

            % decrease maximum value of the energy grid
            eg.updateMaxValue(eg.node(end)/(1+eg.updateFactor));

            % New Solution
            if isnonlinear
                auxEedf = nonLinearSolver(bz);
            else
                auxEedf = bz.linearSolver;
            end

            % check decades of decay
            decades = log10(auxEedf(1)) - log10(auxEedf(end));

        end

    end
    
    % evaluate power balance
    bz.evaluatePower(true);
    
    % evaluate rate coefficients
    bz.evaluateRateCoeff();
    
    % evaluate transport parameters
    bz.evaluateSwarmParameters();
    
    % evaluate first anisotropy
    bz.evaluateFirstAnisotropy();
    
    % bradcast obtention of a solution for the boltzmann equation
    notify(bz, 'obtainedNewEedf');

end