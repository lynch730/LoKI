function obtainTimeDependentSolution(bz)
    % obtainTimeDependentSolutions solves the time-dependent electron
    % Boltzmann equation (to be used for pulsed simulations)
    
    %       % obtain initial solution (steady-state solution for the
    %       initial value of the reduced electric field) if
    %       boltzmann.includeNonConservativeIonization ||
    %       boltzmann.includeNonConservativeAttachment || ...
    %           boltzmann.includeEECollisions
    %         % iterate non-linear solver until convergence
    %         nonLinearSolver = str2func(boltzmann.nonLinearAlgorithm);
    %         boltzmann.isTimeDependent = false;
    %         auxEedf = nonLinearSolver(boltzmann);
    %         boltzmann.isTimeDependent = true;
    %       else
    %         % invert the boltzmann matrix to obtain an eedf (linear
    %         solution) auxEedf = boltzmann.linearSolver;
    %       end
    
    eg = bz.energyGrid;

    % obtain initial solution (as a Maxwellian EEDF at the gas temperature)
    auxEedf = Constant.boltzmannInEV * bz.workCond.gasTemperature;
    auxEedf = exp(-eg.cell ./ auxEedf);
    auxEedf = auxEedf ./ ( eg.step*dot(sqrt(eg.cell), auxEedf) );
    
    % switch for electron density evolution
    bz.eDensIsTimeDependent = false;
    eDen_td_flag = bz.eDensIsTimeDependent;
    
    % save initial value of electron density
    auxNe = bz.workCond.electronDensity;
    
    % evaluate sampling points prescribed by the user
    t0 = bz.pulseFirstStep;
    t1 = bz.pulseFinalTime;
    Nt = bz.pulseSamplingPoints;
    switch bz.pulseSamplingType
        case ('linspace')
            times = linspace(t0, t1, Nt);
        case ('logspace')
            times = logspace(log10(t0), log10(t1), Nt);
    end
    times = [0, times];

    % check if the GUI is activated
    GUIisActive = ~isempty(findobj('name', 'LoKI Simulation Tool'));
    
    % initialise ODE solver parameters

    % flush persistent variables previously stored in function memory
    eedfTimeDerivative(0, [auxEedf auxNe]', bz, true, eDen_td_flag);  

    % load ode options specified by the user
    odeOptionsLocal = bz.odeOptions;
    odeOptionsLocal.NonNegative = 1:length(auxEedf)+1; % ensure >0 EEDF

    % GUI Settings
    if GUIisActive
        odeOptionsLocal.OutputFcn = @odeProgressBar;
        odeProgressBar(0, 0, 'firstInit', t0, t1);
    end
    
    % perform temporal integration along the pulse (dividing integration
    % into smaller chunks for performance sake)
    tdx0 = 1;
    tdx1 = 1;
    finalTimes = [];
    finalEedfs = [];
    finalNe = [];
    
    % Until all time steps hit
    while tdx1 < length(times)

        % evaluate final target time of the current chunk
        if tdx0 == 1
            targetFinalTime = times(2)*10;
        else
            targetFinalTime = times(tdx0)*10;
        end
        
        % find first element of the times array > targetFinalTime
        for idx = tdx1+1:length(times)
            if times(idx) > targetFinalTime
                tdx1 = idx;
                break
            elseif idx == length(times)
                tdx1 = idx;
            end
        end

        % perform temporal integration of the current chunk
        %         [timesSol,eedfsSol,~,~,~] = ...
        [timesSol, solutions] = ode15s(@eedfTimeDerivative, ...
                                     times(tdx0:tdx1), ...
                                     [auxEedf auxNe]', ...
                                     odeOptionsLocal, ...
                                     bz, false, eDen_td_flag);

        % separate solutions
        eedfsSol = solutions(:,1:end-1);
        neSol = solutions(:,end);

        % renormalize eedfs of the current chunck
        for idx = 1:length(timesSol)
            norm = (eg.step*dot(eedfsSol(idx,:), sqrt(eg.cell)));
            eedfsSol(idx,:) = eedfsSol(idx,:)/norm;
        end

        % evaluate initial conditions for the integration of the next chunk
        tdx0 = tdx1;
        auxEedf = eedfsSol(end,:);
        auxNe = neSol(end);
        odeOptionsLocal.InitialStep = timesSol(end)-timesSol(end-1);

        % acumulate solution of the current chuck with the previous ones
        if isempty(finalTimes)
            finalTimes = timesSol;
            finalEedfs = eedfsSol;
            finalNe = neSol;
        else
            finalTimes = [finalTimes; timesSol(2:end)];
            finalEedfs = [finalEedfs; eedfsSol(2:end,:)];
            finalNe = [finalNe; neSol(2:end)];
        end

    end
    
    % close ode integration progress bar in case it was active
    if GUIisActive
        odeProgressBar(0,0,'lastDone');
    end
    
    % evaluate output and save solutions for each time step
    for idx = 1:length(finalTimes)
    
        % evaluate current working conditions(and save them into the
        % working conditions object)
        bz.workCond.currentTime = finalTimes(idx);
        bz.workCond.reducedField = bz.pulseFunction(finalTimes(idx), ...
                                    bz.pulseFunctionParameters);
        bz.workCond.reducedFieldSI = bz.workCond.reducedField * 1.0e-21;
        bz.workCond.electronDensity = finalNe(idx);
    
        % save current eedf in the properties of the Boltzmann object (to
        % be used in other methods)
        bz.eedf = finalEedfs(idx,:);
    
        % evaluate intermidiate quantities needed for power evaluation
        eedfTimeDerivative(finalTimes(idx), ...
                           [finalEedfs(idx,:) finalNe(idx)]', ...
                           bz, false, eDen_td_flag);
    
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

end
