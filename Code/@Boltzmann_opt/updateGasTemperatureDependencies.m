function updateGasTemperatureDependencies(boltzmann, ~, ~)
    % updateGasTemperatureDependencies is a function that evaluates all the
    % dependencies of the boltzmann object (either direct or indirect) on
    % the gasTemperature property of the workCond object
    
    % evaluate population of state dependent on the gas temperature (e.g.
    % Boltzmann distributions)
    populationsDependentOnGasTemperature = false;
    for gas = boltzmann.gasArray
        for state = gas.stateArray
            if ~isempty(state.populationFunc)
                for parameter = state.populationParams
                    if strcmp(parameter{1}, 'gasTemperature')
                        state.evaluatePupulation(boltzmann.workCond);
                        populationsDependentOnGasTemperature = true;
                        break;
                    end
                end
            end
        end
    end
    
    % evaluate Boltzmann operators dependent on the gas temperature
    if populationsDependentOnGasTemperature 
        boltzmann.evaluateMatrix;
    else
        boltzmann.evaluateContinuousOperators;
    end
    
end