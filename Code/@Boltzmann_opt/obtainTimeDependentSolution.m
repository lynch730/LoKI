function obtainTimeDependentSolution(bz)
    % obtainTimeDependentSolutions solves the time-dependent electron
    % Boltzmann equation (to be used for pulsed simulations)
    
    eg = bz.energyGrid;
    N = eg.cellNumber;

    %% Time Settings

    % switch for electron density evolution
    bz.eDensIsTimeDependent = false;
    eDen_td_flag = bz.eDensIsTimeDependent;
    
    % evaluate sampling points prescribed by the user
    t1 = bz.pulseFirstStep;
    t2 = bz.pulseFinalTime;
    Nt = bz.pulseSamplingPoints;
    switch bz.pulseSamplingType
        case ('linspace')
            times = linspace(t1, t2, Nt);
        case ('logspace')
            times = logspace(log10(t1), log10(t2), Nt);
    end
    times = [0, times];
    Nt = numel(times);

    % check if the GUI is activated
    GUIisActive = ~isempty(findobj('name', 'LoKI Simulation Tool'));
    

    %% Initialise ODE solver parameters

    % obtain initial solution (as a Maxwellian EEDF at the gas temperature)
    eedf = Constant.boltzmannInEV * bz.workCond.gasTemperature;
    eedf = exp(-eg.cell ./ eedf);
    eedf = eedf ./ ( eg.step*dot(sqrt(eg.cell), eedf) );

    % save initial value of electron density
    ne = bz.workCond.electronDensity;
    
    % Store Parameters in Memory
    p = init_td_param(bz, eDen_td_flag);
    
    % Split times on log-scale by exponent
    [~, tdx] = unique(ceil(log10(times(2:end))));
    tdx = [1; tdx+1; Nt]; % Include first and last time steps
    tdx(logical([0; diff(tdx)==1])) = []; % Clear short steps (or will mess up return)
    Ngroup = numel(tdx)-1; 

    % Set ode function with parameters p
    ode_func = @(t, y) eedfTimeDerivative_opt(t, y, p, 0);
    jac_func = @(t, y) eedfTimeDerivative_opt(t, y, p, 1);

    % load ode options specified by the user
    odeOptionsLocal = bz.odeOptions;
%     odeOptionsLocal.NonNegative = 1:length(eedf)+1; % ensure >0 EEDF
    odeOptionsLocal.Jacobian = jac_func;
    odeOptionsLocal.AbsTol = 1.0e-10;
    odeOptionsLocal.RelTol = 1.0e-5;
    odeOptionsLocal.NormControl = 'off';
%     odeOptionsLocal.BDF = 'on';
%     odeOptionsLocal.MaxOrder = 5;
    odeOptionsLocal.Stats = 'on';

    % GUI Settings
    if GUIisActive
        odeOptionsLocal.OutputFcn = @odeProgressBar;
        odeProgressBar(0, 0, 'firstInit', t1, t2);
    end

    % Insert first
    sol = zeros(Nt, N+1);
    sol(1,:) = [eedf, ne];

    %% Main Loop
    for i = 1:Ngroup
        
        % Index of Current Chunk
%         time_ind = tdx(i)-(i>1):tdx(i+1);
        time_ind = tdx(i):tdx(i+1);

        % Run ODE
        [~, sol_tmp] = ode15s(ode_func, ...
                         times(time_ind), ...
                         sol(tdx(i), :)', ...
                         odeOptionsLocal);
        
        % separate solutions
        sol(time_ind, :) = sol_tmp;

        % renormalize eedfs of the current chunck
        norm = eg.step*( sol_tmp(:, 1:N) * sqrt(eg.cell(:)) );
        sol(time_ind, 1:N) = sol(time_ind, 1:N) ./ norm;

    end
    
    % close ode integration progress bar in case it was active
    if GUIisActive
        odeProgressBar(0,0,'lastDone');
    end
    
    % Temp plot
    figure(1); clf; hold on;
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    ax.XLim = [eg.cell(1), eg.cell(end)];
    ax.YLim = [1e-100 1e20];
    pp = plot(eg.cell, sol(1, 1:N), '-k');
    for i = 2:Nt
        pp.YData = abs(sol(i, 1:N));
        drawnow
        pause(0.01)
    end
    
    % evaluate output and save solutions for each time step
    for i = 1:Nt
        
        % evaluate current working conditions(and save them into the
        % working conditions object)
        bz.workCond.currentTime = times(i);
        bz.workCond.reducedField = bz.pulseFunction(times(i), ...
                                    bz.pulseFunctionParameters);
        bz.workCond.reducedFieldSI = bz.workCond.reducedField * 1.0e-21;
        bz.workCond.electronDensity = sol(i, N+1);
        
        % save current eedf in the properties of the Boltzmann object (to
        % be used in other methods)
        bz.eedf = sol(i, 1:N);
        
        % evaluate intermidiate quantities needed for power evaluation
        eedfTimeDerivative(times(i),...
                           [sol(i, 1:N) sol(i, N+1)]', ...
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


function p = init_td_param(bz, eDen_td_flag)
    
    % Store Constants
    p.eDen_td_flag = eDen_td_flag;
    p.includeNonConservativeIonization = bz.includeNonConservativeIonization;
    p.includeNonConservativeAttachment = bz.includeNonConservativeAttachment;
    p.eDensGrowthModel = bz.eDensGrowthModel;
    p.includeEECollisions = bz.includeEECollisions;
    p.ne_gasDensity = bz.workCond.gasDensity;

    con = Constant();
    p.qe = con.electronCharge;                % electron charge (SI units)
    p.gamma = con.gamma;                     % p.gamma parameter (sqrt(2e/me))
    p.eps0 = con.vacuumPermittivity;           % vacuum permitivity (SI units)
    
    % Store Energy Grid Data
    p.N = bz.energyGrid.cellNumber;   % number of cells used in the energy grid
    p.deps = bz.energyGrid.step;         % energy step of the energy grid
    p.eCell = bz.energyGrid.cell;         % energy at cells of the energy grid
    p.eNode = bz.energyGrid.node;         % energy at nodes of the energy grid
    p.Nden = bz.workCond.gasDensity;              % gas density (SI units)
    p.WoN = bz.workCond.reducedExcFreqSI;      % reduced angular exitation frequency (SI units)
    
    % Store E-Field Time Dependence
    if bz.isTimeDependent
        p.EoN = @(t) bz.pulseFunction(t, bz.pulseFunctionParameters).*1e-21;
    else
        p.EoN = @(t) bz.workCond.reducedFieldSI + 0.0.*t;
    end

    % Store Matrices

    % evaluate basic boltzmann matrix (without field, ionization,
    % attachment or e-e collisions operators)
    p.M = bz.elasticMatrix + bz.CARMatrix + bz.inelasticMatrix;
    
    % save local copy of the regular field operator (without growth models)
    p.Mfield = bz.fieldMatrix;

    % include ionization operator (either conservative or not)
    p.Mncon = zeros(p.N);
    if bz.includeNonConservativeIonization
        p.M = p.M + bz.ionizationMatrix;
        p.Mncon = p.Mncon + bz.ionizationMatrix;
    else
        p.M = p.M + bz.ionizationConservativeMatrix;
    end

    % include attachment operator (either conservative or not)
    if bz.includeNonConservativeAttachment
        p.M = p.M + bz.attachmentMatrix;
        p.Mncon = p.Mncon + bz.attachmentMatrix;
    else
        p.M = p.M + bz.attachmentConservativeMatrix;
    end

    % evaluate time independent elements of the growth model operators
    % (in case they are activated)
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        
        % evaluate integrand to calculate the effective ionization rate
        p.CIEffIntegrand = sum(p.Mncon)';

        % local copies of total momentum transfer cross sections (at
        % nodes and cells)
        p.sig_total = bz.totalCrossSection;
        p.cTCS = 0.5*(p.sig_total(1:end-1) + p.sig_total(2:end));
        
        % choose growth model for the electron density (either spatial
        % or temporal)
        switch bz.eDensGrowthModel
            case 'temporal'

                % evaluate M elements of the temporal growth
                % operator (e-density time variation term, dn/dt)
                p.M_Gdiag = -sqrt(p.eCell)/p.gamma;

            case 'spatial'

                % evaluate the diffusion and mobility components of the
                % spatial growth operator
                p.M_Gdiag = p.eCell./(3*p.cTCS);

                p.M_GMSE = ...     % mobility component (diag sup)
                    1/(6*p.deps)*[0 p.eCell(1:p.N-1)./(p.cTCS(1:p.N-1))];
                p.M_GMINFE = ...     % mobility component (diag inf)
                    -1/(6*p.deps)*[p.eCell(2:p.N)./(p.cTCS(2:p.N)) 0];
                
                % evaluate components of the extra electric field operator of the spatial growth model
                p.M_eFSG = p.eNode/p.deps./(6*p.sig_total);
                p.M_eFSG(1) = 0;
                p.M_eFSG(end) = 0;
            
        end

    end

    % evaluate time idependent elements of the e-e collision operator
    % (in case that e-e collisions are activated)
    if bz.includeEECollisions
        
        % evaluating auxiliary matrix used in the evaluation of the ee
        % upflux vector without multiplicative constant
        p.eeMatrixAuxA = zeros(p.N);
        auxEnergyArray = -(p.deps/2)*sqrt(p.eCell)+(2/3)*p.eCell.^1.5;
        
        % because of power conservation, terms on the last row
        % (p.N,:) and on the first column (:,1) are zero
        for k=1:p.N-1
            p.eeMatrixAuxA(k,2:k) = auxEnergyArray(2:k);
            p.eeMatrixAuxA(k,k+1:p.N) = (2/3)*p.eNode(k+1)^1.5;
        end
        
        % detailed balance condition
        for k=1:p.N-1
            p.eeMatrixAuxA(k,2:p.N) = sqrt(p.eeMatrixAuxA(k,2:p.N).*p.eeMatrixAuxA(1:p.N-1,k+1)');
        end
        
        % evaluating auxiliary matrix used in the evaluation of the ee downflux vector without multiplicative constant
        p.eeMatrixAuxB = transpose(p.eeMatrixAuxA);
        
    end
    
    % Convert to sparse 
    p.M = sparse(p.M) ./ sqrt(p.eCell');

end
