function eedf = mixingDirectSolutions(bz)
    % mixingDirectSolutions solves the non-linear Boltzmann equation (i.e.
    % including one of more of the following operators: non-conservative
    % ionization, non-conservative attachment or e-e collision) with an
    % algorithm of mixing direct solutions of a linearized Boltzmann
    % equation.
    
    % local copies of energy grid variables
    N = bz.energyGrid.cellNumber;
    
    % initial guess for the eedf from the linear Boltzmann equation (i.e.
    % without non-linear operators)
    bz.linearSolver();
    
    % select the appropiate method depending on the growth model (if
    % non-conservative operators are activated)
    solveEGrowthModel = [];
    if bz.includeNonConservativeIonization || ...
        bz.includeNonConservativeAttachment
        switch bz.eDensGrowthModel
            case 'spatial'
                solveEGrowthModel = str2func('solveSpatialGrowthMatrix');
                % reset spatial growth related matrices
                bz.ionSpatialGrowthD  = zeros(N);
                bz.ionSpatialGrowthU = zeros(N);
                bz.fieldMatrixSpatGrowth = zeros(N);
            case 'temporal'
                solveEGrowthModel = str2func('solveTemporalGrowthMatrix');
                % reset temporal growth related matrices
                bz.fieldMatrixTempGrowth = zeros(N);
                bz.ionTemporalGrowth = zeros(N);
        end
    end
    
    % solve the non-linear Boltzmann system to obtain an eedf (different
    % cases depending on the activated operators)
    if bz.includeEECollisions

        % reset e-e collisions related variables
        bz.alphaee = 0;
        bz.Aee = zeros(N);
        bz.Bee = zeros(N);
        bz.Afinal = zeros(1, N);
        bz.Bfinal = zeros(1, N);

        % Switch growth model
        if isempty(solveEGrowthModel)
            eedf = bz.solveEEColl(); % solution with the e-e collisions only
        else

            % counter for the number of iterations of the global cycle
            globalIteration = 0;        

            % maximum number of iterations for the global cycle
            maxGlobalIterations = 20;  
            
            % solution with e-e collisions and growth model
            while globalIteration < maxGlobalIterations

                % cycle updating only the growth model and EE operator
                eedfOld = solveEGrowthModel(bz);   
                eedf = bz.solveEEColl();

                % exit global cycle in case of convergence (maximum
                % relative difference of each value of the eedf < 1e-9)
                if max(abs(eedf-eedfOld)./eedfOld) < bz.maxEedfRelError
                    break;
                end

                % Advance Iteraiton
                globalIteration = globalIteration+1;

                % throw a warning message if the global cycle has not converged
                if globalIteration == maxGlobalIterations
                    warning(['Global cycle (mixing of direct solutions): ',...
                            'EEDF did not converge after %d iterations\n'], ...
                            maxGlobalIterations);
                end

            end
        end
        
    else
        % solution with growth model only
        eedf = solveEGrowthModel(bz);
    end

end