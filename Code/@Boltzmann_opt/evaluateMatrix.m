function evaluateMatrix(boltzmann, ~, ~)
    
    % evaluate continuous operators of the Boltzmann equation
    % (field+elastic+CAR operators)
    boltzmann.evaluateContinuousOperators();
    
    % evaluate discrete operators of the Boltzmann equation (inelastic and
    % superelastic operators without including the contributions from
    % ionization nor attachment)
    boltzmann.evaluateInelasticOperators();
    
    % evaluate ionization operator of the Boltzmann equation
    boltzmann.evaluateIonizationOperator();
    
    % evaluate attachment operator of the Boltzmann equation
    boltzmann.evaluateAttachmentOperator();

end