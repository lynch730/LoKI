function evaluateContinuousOperators(boltzmann)
    
    % evaluate total and elastic momentum transfer cross sections for
    % continuous operators
    boltzmann.evaluateTotalAndElasticCrossSections();
    
    % evaluate field operator
    boltzmann.evaluateFieldOperator();
    
    % evaluate elastic operator
    boltzmann.evaluateElasticOperator();
    
    % evaluate CAR operator
    if ~isempty(boltzmann.CARgases)
        boltzmann.evaluateCAROperator();
    end

end