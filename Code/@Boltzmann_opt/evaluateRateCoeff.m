function [rateCoeffAll, rateCoeffExtra] = evaluateRateCoeff(bz)  
    % initialize rateCoeffAll and rateCoeffExtra structures
    rateCoeffAll = struct.empty;
    rateCoeffExtra = struct.empty;
    
    % evaluate rate coefficient for all collision
    for gas = bz.gasArray
        % collisions taken into account for solving the eedf
        for collision = gas.collisionArray
            rateCoeffAll(end+1).collID = collision.ID;
            [ineRate, supRate] = collision.evaluateRateCoeff(bz.eedf);
            rateCoeffAll(end).value = [ineRate, supRate];
            rateCoeffAll(end).energy = collision.threshold;
            rateCoeffAll(end).collDescription = collision.description;
        end
        % collisions not taken into account for solving the eedf
        for collision = gas.collisionArrayExtra
            rateCoeffExtra(end+1).collID = collision.ID;
            [ineRate, supRate] = collision.evaluateRateCoeff(bz.eedf);
            rateCoeffExtra(end).value = [ineRate, supRate];
            rateCoeffExtra(end).energy = collision.threshold;
            rateCoeffExtra(end).collDescription = collision.description;
        end
    end
    
    % store rate coefficients information in the boltzmann properties
    bz.rateCoeffAll = rateCoeffAll;
    bz.rateCoeffExtra = rateCoeffExtra;

end