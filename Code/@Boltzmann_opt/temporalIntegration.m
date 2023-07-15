function eedf = temporalIntegration(bz)
    % eedfTimeIntegration performs a time integration of the Boltzmann
    % equation (to be completed when the function is completely developed)
    
    % initial guess for the eedf from the linear Boltzmann equation (i.e.
    % without non-linear operators)
    eedf = bz.linearSolver();
    
    % iterative ode solver
    ne0 = bz.workCond.electronDensity;
    eedfTimeDerivative(0, [eedf ne0]', bz, true, false); % flush persistent

    % Fill ODE Options
    odeOptionsLocal = bz.odeOptions;
    odeOptionsLocal.NonNegative = 1:numel(eedf);
    odeOptionsLocal.Events = @steadyStateEventFcn;
    
    % Run ODE15s
    [~,~,~,variablesSolution,~] = ode15s(@eedfTimeDerivative, ...
                                          [0 inf], ...
                                          [eedf ne0]', ...
                                          odeOptionsLocal, ...
                                          bz, 0, 0);
    
    % select final solution if more than one is found
    if length(variablesSolution(:,1))>1
        eedf = variablesSolution(end, 1:end-1);
    else
        eedf = variablesSolution(1, 1:end-1);
    end
    
    % ensure proper normalization
    norm = sum(eedf.*sqrt(bz.energyGrid.cell)) * bz.energyGrid.step;
    eedf = eedf/norm;
    
    % store eedf in the properties of the boltzmann object
    bz.eedf = eedf;

end
