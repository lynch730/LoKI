function evaluateMacroscopicParameters(bz)
    % evaluateMacroscopicParameters evaluate elctron-impact rate
    %   coefficients, swarm parameters and power channels using Boltzmann
    %   operators from the Boltzmann equation
    
    bz.evaluatePower(false);
    bz.evaluateSwarmParameters;
    bz.evaluateRateCoeff;

end