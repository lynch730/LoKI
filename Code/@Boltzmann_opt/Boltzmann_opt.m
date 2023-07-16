classdef Boltzmann_opt < handle
    %Boltzmann Class that solves the boltzmann equation under certain conditions to
    %   obtain an EEDF
    %
    %   Boltzmann creates the different boltzmann matrices corresponding
    %   to the different elements of the boltzmann equation. In order to
    %   create a BoltzmannEq object, a gasArray with the gas mixture
    %   information and an energyGrid must be sent to the constructor.
    
    %% Property Definition
    properties
        
        gasArray = Gas.empty;                   % handle array to the gas mixture
        energyGrid = Grid.empty;                % handle to the energy grid in which the Boltzmann equation is to be solved
        workCond = WorkingConditions.empty;     % handle to the working conditions of the simulation
        CARgases = [];                          % gases described by the continuous approximation for rotations (CAR)
        
        ionCollOpType;            % describes the type of ionization collision operator
        eDensGrowthModel;         % describes the electron density growth model
        includeEECollisions;      % indicates if electron-electrons collisions are accounted for
        
        nonLinearAlgorithm;       % algorithm used to solve the non-linear operators (ionization, attachment or e-e)
        mixingParameter;          % solution mixing fraction used on the iterative scheme of the growth routines
        maxEedfRelError;          % maximum relative difference for the eedf between two consecutive iterations (non-linear routines)
        maxPowerBalanceRelError;  % maximum value for the relative power balance (threshold for the warning message)
        odeOptions;               % user defined options for the ODE solver (used when nonLinearAlgorithm is temporalIntegration)
    
        isTimeDependent = false;        % false when solving quasi-stationary Boltzmann equation, true for time-dependent pulsed Boltzmann equation
        pulseFunction;                  % function handle to the reduced electric field pulse function
        pulseFirstStep = [];            % first time step for the pulse (only for the sampling of the solution)
        pulseFinalTime = [];            % final time for the pulse
        pulseSamplingType = [];         % type of sampling for time (either 'linspace' or 'logspace')
        pulseSamplingPoints = [];       % number of sampling points for the tinal solution
        pulseFunctionParameters = {};   % additional parameters to be sent to the reduced electric field pulse function
        eDensIsTimeDependent = false;   % true if electron density should evolve in pulsed simulations
        
        totalCrossSection = [];   % total momentum transfer cross section
        elasticCrossSection = []; % total elastic cross section
        
        g_E = [];                   % elements of the field operator
        g_c = [];                   % elements of the elastic operator
        g_car = [];                 % elements of the CAR operator
        g_fieldSpatialGrowth = [];  % elements of the field operator with spatial electron density growth model
        g_fieldTemporalGrowth = []; % elements of the field operator with temporal electron density growth model
        Afinal = [];                % elements of the energy upflux from electron-electron collisions (used in power balance)
        Bfinal = [];                % elements of the energy downflux from electron-electron collisions (used in power balance)
        
        ionizationThresholdIsSmallerThanUmax;     % indicates if at least one ionization threshold is within the energy grid range
        includeNonConservativeIonization = false; % boolean that indicates if growth models are to be used because of ionization
        includeNonConservativeAttachment = false; % boolean that indicates if growth models are to be used because of attachment
        attachmentThresholdIsSmallerThanUmax;     % indicates if at least one attachment threshold is within the energy grid range
        
        CIEff;                    % effective ionization rate coefficient
        alphaRedEff;              % reduced first Townsend ionization coefficient
        alphaee;                  % alpha constant of electron-electron collisions
        
        fieldMatrix = [];           % matrix of the electric field operator of the boltzmann equation (continuous term)
        fieldMatrixTempGrowth = []; % matrix of the electric field operator of the boltzmann equation with temporal growth model
        fieldMatrixSpatGrowth = []; % matrix of another part of the electric field operator of the boltzmann equation with spatial growth model
        elasticMatrix = [];         % matrix of the elastic collision operator of the boltzmann equation (continuous term)
        CARMatrix = [];             % matrix of the CAR operator of the boltzmann equation (continuous term)
        continuousMatrix = [];      % matrix of the continuous operator of the boltzmann equation (sum of all continuos terms)
        inelasticMatrix = [];       % matrix of the discrete operator of the boltzmann equation (sum of all discrete terms)
        ionizationConservativeMatrix = []; % matrix of the conservative ionization operator of the boltzmann equation (discrete terms)
        ionizationMatrix = [];      % matrix of the non-conservative ionization operator of the boltzmann equation (discrete terms)
        attachmentConservativeMatrix = []; % matrix of the conservative attachment operator of the boltzmann equation (discrete terms)
        attachmentMatrix = [];      % matrix of the non-conservative attachment operator of the boltzmann equation (discrete terms)
        ionTemporalGrowth= [];      % matrix of the electron-density temporal growth term
        ionSpatialGrowthD = [];     % matrix of the electron-density spatial growth diffusion term
        ionSpatialGrowthU = [];     % matrix of the electron-density spatial growth mobility term
        Aee = [];                   % auxiliary matrix to calculate upflux vector of electron-electron collisions
        Bee = [];                   % auxiliary matrix to calculate downflux vector of electron-electron collisions
        
        eedf = [];                    % eedf obtained as solution to the boltzmann equation
        power = struct.empty;         % power balance of the boltzmann equation
        swarmParam = struct.empty;    % swarm parameters obtained with the eedf
        rateCoeffAll = struct.empty;  % rate coefficients obtained with the eedf and collisions in gasArray
        rateCoeffExtra = struct.empty;% extra rate coefficients for collisions not taken into account to obtain the eedf
        firstAnisotropy = [];         % value of the first anisotropy of the EDF (two term approximation)
        
    end

    %% events
    events
        genericStatusMessage;
        obtainedNewEedf;
    end

    %% Public Methods
    methods (Access = public)
        
        %% Constructor
        function bz = Boltzmann_opt(setup)
                    
            % store the gas array
            bz.gasArray = setup.electronKineticsGasArray;
            bz.energyGrid = setup.energyGrid;
            
            % store the energy grid and add corresponding listener
            addlistener(bz.energyGrid, 'updatedMaxEnergy2', @bz.evaluateMatrix);
            
            % store working conditions and add corresponding listeners
            bz.workCond = setup.workCond;
            addlistener(bz.workCond, 'updatedGasTemperature', @bz.evaluateMatrix);
            %       addlistener(workCond, 'updatedGasDensity', @boltzmann.);
            addlistener(bz.workCond, 'updatedReducedField', @bz.evaluateFieldOperator);
            addlistener(bz.workCond, 'updatedExcitationFrequency', @bz.evaluateFieldOperator);
            
            % store gases for which the CAR is activated (in case there is any)
            ek = setup.info.electronKinetics;
            if isfield(ek, 'CARgases')
                bz.CARgases = ek.CARgases;
            end
            
            % store electron-impact ionization operator type
            bz.ionCollOpType = ek.ionizationOperatorType;
            
            % store electron density growth model type
            bz.eDensGrowthModel = ek.growthModelType;
            
            % store electron-electron collisions setups
            bz.includeEECollisions = ek.includeEECollisions;
            
            % store information about the numerical details on how to solve
            % the Boltzmann equation
            ek_nlr = ek.numerics.nonLinearRoutines;
            bz.nonLinearAlgorithm = ek_nlr.algorithm;
            bz.maxEedfRelError = ek_nlr.maxEedfRelError;
            if strcmp(bz.nonLinearAlgorithm, 'mixingDirectSolutions')
                bz.mixingParameter = ek_nlr.mixingParameter;
            else
                % store configuration of the ODE solver
                options = odeset();
                for parameter = fields(options)'
                    if isfield(ek_nlr, 'odeSetParameters') && ...
                            isfield(ek_nlr.odeSetParameters, parameter{1})
                        options.(parameter{1}) = ek_nlr.odeSetParameters.(parameter{1});
                    end
                end
                bz.odeOptions = options;
            end
            bz.maxPowerBalanceRelError = ek.numerics.maxPowerBalanceRelError;
            
            % check if the the simulation is steady-state or pulsed
            if setup.pulsedSimulation
                
                % store information about pulsed simulation if activated
                bz.isTimeDependent = true;
                pinfo = setup.pulseInfo;
                bz.pulseFunction = pinfo.function;
                bz.pulseFirstStep = pinfo.firstStep;
                bz.pulseFinalTime = pinfo.finalTime;
                bz.pulseSamplingType = pinfo.samplingType;
                bz.pulseSamplingPoints = pinfo.samplingPoints;
                bz.pulseFunctionParameters = pinfo.functionParameters;
            
                % set initial value of the reduced electric field in the
                % working conditions object
                bz.workCond.reducedField = bz.pulseFunction(0, bz.pulseFunctionParameters);
                bz.workCond.reducedFieldSI = bz.workCond.reducedField*1e-21;
            
            end
            
            % % allocate memory for different properties
            % boltzmann.totalCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
            % boltzmann.elasticCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
            % boltzmann.fieldMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.elasticMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            bz.CARMatrix = zeros(bz.energyGrid.cellNumber,bz.energyGrid.cellNumber);
            % boltzmann.g_car = zeros(1,boltzmann.energyGrid.cellNumber+1);
            % boltzmann.continuousMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.discreteMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.ionizationIneMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            
            % % initialize quantities of the ionization and electron-electron routines
            % boltzmann.ionizationMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.Aee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.Bee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
            % boltzmann.ionSpatialGrowthD  = zeros(boltzmann.energyGrid.cellNumber);
            % boltzmann.ionSpatialGrowthU = zeros(boltzmann.energyGrid.cellNumber);
            % boltzmann.fieldMatrixSpatGrowth = zeros(boltzmann.energyGrid.cellNumber);
            % boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber);
            % boltzmann.ionTemporalGrowth = zeros(boltzmann.energyGrid.cellNumber);

            bz.evaluateMatrix(); % Evaluate boltzmann matrix

        end
        
        % Declerations
        solve(boltzmann)
        evaluateMacroscopicParameters(boltzmann)
        updateDensityDependencies(boltzmann)
        
    end
    
    %% Private Methods
    methods (Access = private)
        
        % Helper
        bz = process_setup(setup)
        updateGasTemperatureDependencies(boltzmann, ~, ~)
        
        % Operators
        evaluateMatrix(boltzmann, ~, ~)
        evaluateContinuousOperators(boltzmann)
        evaluateTotalAndElasticCrossSections(boltzmann)
        evaluateFieldOperator(boltzmann, ~, ~)
        evaluateElasticOperator(boltzmann)
        evaluateCAROperator(boltzmann)
        evaluateInelasticOperators(boltzmann)
        evaluateIonizationOperator(boltzmann)
        evaluateAttachmentOperator(boltzmann)
         
        % Wrappers
        obtainTimeIndependentSolution(boltzmann)
        obtainTimeDependentSolution(boltzmann)
    
        % EEDF Solvers
        eedf = linearSolver(boltzmann)
        eedf = temporalIntegration(boltzmann)
        eedf = mixingDirectSolutions(boltzmann)
        eedf = solveSpatialGrowthMatrix(boltzmann)
        eedf = solveTemporalGrowthMatrix(boltzmann)
        eedf = solveEEColl(boltzmann)
          
        % Post Processing
        power = evaluatePower(boltzmann, checkPowerBalance)
        swarmParam = evaluateSwarmParameters(boltzmann)
        [rateCoeffAll, rateCoeffExtra] = evaluateRateCoeff(boltzmann)
        evaluateFirstAnisotropy(boltzmann)
    
    end
  
end

