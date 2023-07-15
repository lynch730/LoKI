function eedf = solveEEColl(bz)
    
    e = Constant.electronCharge;                % electron charge
    e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
    ne = bz.workCond.electronDensity;    % electron density (SI units)
    n0 = bz.workCond.gasDensity;         % gas density (SI units)
    EoN = bz.workCond.reducedFieldSI;    % reduced electric field (SI units)
    EoNTd = bz.workCond.reducedField;    % reduced electric field (Td)

    energyCell = bz.energyGrid.cell;
    energyStep = bz.energyGrid.step;
    N = bz.energyGrid.cellNumber;
    energyNode = bz.energyGrid.node;
    eedf = bz.eedf;
    
    Mee = zeros(N);
    
    % writting of Boltzmann matrix without the electron-electron collisional operator
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        switch bz.eDensGrowthModel
            case 'spatial'
                baseMatrix = bz.ionizationMatrix + bz.attachmentMatrix + bz.elasticMatrix + ...
                    bz.CARMatrix + bz.inelasticMatrix + EoN^2*bz.fieldMatrix + ...
                    bz.ionSpatialGrowthD + bz.ionSpatialGrowthU +  bz.fieldMatrixSpatGrowth;
            case 'temporal'
                baseMatrix = bz.ionizationMatrix + bz.attachmentMatrix + bz.elasticMatrix + ...
                    bz.CARMatrix + bz.inelasticMatrix + bz.fieldMatrixTempGrowth + ...
                    bz.ionTemporalGrowth;
        end
    else
        baseMatrix = bz.ionizationConservativeMatrix + bz.attachmentConservativeMatrix + ...
            bz.elasticMatrix + bz.CARMatrix + bz.inelasticMatrix + EoN^2*bz.fieldMatrix ;
    end
    
    % writing auxiliary matrix A, without constant alpha
    auxA = zeros(N);
    auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
    % from power conservation, terms on the last row (cellNumber,:) and on the first column (:,1) are zero
    for k=1:N-1
        auxA(k,2:k) = auxEnergyArray(2:k);
        auxA(k,k+1:N) = (2/3)*energyNode(k+1)^(3/2);
    end
    % detailed balance condition
    for k=1:N-1
        auxA(k,2:N) = sqrt(auxA(k,2:N).*auxA(1:N-1,k+1)');
    end
    
    % writing auxiliary matrix B, without constant alpha
    auxB = transpose(auxA);
    
    % initial value for eedfOld (calculated without electron-electron collisional operator)
    eedfOld = eedf;
    % initialize cycle counter
    iteration = 0;
    
    while true
        % Calculation of constant alpha
        meanEnergy = energyStep*dot(energyCell.^(3/2),eedf);  % mean energy
        Te = (2/3)*meanEnergy;                                % electron temperature in eV Te = (2/3)*meanEnergy
        logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));           % Coulomb logarithm
        alpha = (ne/n0)*(e^2/(8*pi*e0^2))*logC;               % alpha
    
        % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
        A = (alpha/energyStep)*(auxA*eedf');
        B = (alpha/energyStep)*(auxB*eedf');
    
        % writing of the electron-electron operator
        Mee(1:N+1:N*N) = -(A(1:N)+B(1:N));
        Mee(N+1:N+1:N*N) = B(2:N);
        Mee(2:N+1:N*N) = A(1:N-1);
    
        % sum all contributions to the Boltzmann matrix and reescale it to avoid very small numbers
        matrixAux = 1e20*(baseMatrix + Mee);
    
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(matrixAux, bz.energyGrid);
    
        % calculation of ratio = Pee/PRef
        bz.eedf = eedf;
        bz.Afinal = A;
        bz.Bfinal = B;
        bz.evaluatePower(false);
        Preference = bz.power.reference;
        Pee = bz.power.electronElectronNet;
        ratio = abs(Pee/Preference);
    
        % evaluate convergence criteria
        if  max(abs(eedf-eedfOld)./eedfOld)<bz.maxEedfRelError && abs(ratio)<1e-9
            break;
        elseif max(abs(eedf-eedfOld)./eedfOld)<bz.maxEedfRelError && iteration>200
            str = sprintf('\\t- e-e iterative scheme: EEDF has converged but abs(Pee/Pref)= %.16g > 1e-9\\n',ratio);
            notify(bz, 'genericStatusMessage', StatusEventData(str, 'status'));
            str = sprintf('\\t- Ionization degree:%f; Reduced electric field:%f (Td) \\n',ne/n0,EoNTd);
            notify(bz, 'genericStatusMessage', StatusEventData(str, 'status'));
            break;
        elseif iteration == 400
            warning('Electron-electron iterative scheme: EEDF did not converge\n');
            break;
        end
    
        % update "old" values
        if iteration == 0
            eedfOld = eedf;
        else  % acceleration scheme based on the Newton-Raphson
            % eedf estimate (avoid negative values)
            eedfAux = abs(eedf - (eedf-eedfOld)*ratio/(ratio-ratioOld));
            eedfOld = eedf;
            eedf = eedfAux;
        end
        ratioOld = ratio;
        % update cycle counter
        iteration = iteration+1;
    
    end
    
    % saving auxiliary matrix to be used on the growth model routine
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        bz.Aee = auxA;
        bz.Bee = auxB;
        bz.alphaee = alpha;
    end
    
    % saving eedf
    bz.eedf = eedf;

end