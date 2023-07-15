function [value,isTerminal,direction] = steadyStateEventFcn(~,eedf,boltzmann,~,~)
    % evaluate time derivative of each discrete component of the eedf
    
    persistent eedfOld;
    persistent maxEedfDifference;
    if isempty(eedfOld)
        eedfOld = eedf;
        maxEedfDifference = boltzmann.maxEedfRelError;
        value = 1;
    elseif max(abs(eedf-eedfOld)./eedfOld) < maxEedfDifference
        value = 0;
    else
        eedfOld = eedf;
        value = 1;
    end
    
    isTerminal = 1;   % stop the integration
    direction = 0;    % stop in any direction
    
end