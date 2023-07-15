function power = evaluatePower(bz, checkPowerBalance)
    
    % initialize power structure
    power = struct('field', 0, 'elasticNet', 0, 'elasticGain', 0, 'elasticLoss', 0, 'carNet', 0, 'carGain', 0, ...
        'carLoss', 0, 'excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
        'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
        'ionizationIne', 0, 'attachmentIne', 0, 'inelastic', 0, 'superelastic', 0, 'eDensGrowth', 0, ...
        'electronElectronNet', 0, 'electronElectronGain', 0, 'electronElectronLoss', 0, 'gases', '');
    
    % save a local copy of the EEDF because of performance reasons
    eedfLocal = bz.eedf;
    
    % save a local copy of the energy grid information
    energyGridLocal = bz.energyGrid;   % energyGrid object
    N = energyGridLocal.cellNumber;           % number of cells in the energy grid
    energyStep = energyGridLocal.step;        % energy step
    energyCell = energyGridLocal.cell;        % value of the energy at the cells
    energyNode = energyGridLocal.node;        % value of the energy at the nodes
    
    % multiplicative constant to obtain the right units
    gamma = Constant.gamma;
    
    % auxiliary quantities needed to evaluate the elastic and CAR powers
    kTg = Constant.boltzmannInEV*bz.workCond.gasTemperature;
    aux1 = kTg+energyStep*0.5;
    aux2 = kTg-energyStep*0.5;
    
    % evaluate power absorved per electron at unit gas density due to elastic collisions
    g_cLocal = bz.g_c; % elements of the elastic collision operator (local copy)
    power.elasticNet = gamma*sum(eedfLocal.*(g_cLocal(2:end)*aux2-g_cLocal(1:end-1)*aux1));
    power.elasticGain = gamma*kTg*sum(eedfLocal.*(g_cLocal(2:end)-g_cLocal(1:end-1)));
    power.elasticLoss = power.elasticNet-power.elasticGain;
    
    % evaluate power absorved per electron at unit gas density due to rotations CAR
    if ~isempty(bz.CARgases)
        g_carLocal = bz.g_car; % elements of the CAR operator (local copy)
        power.carNet = gamma*sum(eedfLocal.*(g_carLocal(2:end)*aux2-g_carLocal(1:end-1)*aux1));
        power.carGain = gamma*kTg*sum(eedfLocal.*(g_carLocal(2:end)-g_carLocal(1:end-1)));
        power.carLoss = power.carNet-power.carGain;
    end
    
    % evaluate power gained from electric field and lost due to electron density growth terms
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        switch bz.eDensGrowthModel
            case 'temporal'
                g_ELocal = bz.workCond.reducedFieldSI^2*bz.g_fieldTemporalGrowth; % elements of the electric field operator (local copy)
                power.field = gamma*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
    
                power.eDensGrowth = - bz.CIEff*energyStep*sum(eedfLocal.*energyCell.*sqrt(energyCell));
            case 'spatial'
                totalCrossSectionLocal = bz.totalCrossSection;
                cellTotalCrossSection =  0.5*(totalCrossSectionLocal(1:N) + totalCrossSectionLocal(2:N+1));
                alphaRedEffLocal = bz.alphaRedEff;
                reducedFieldSILocal = bz.workCond.reducedFieldSI;
                % elements of the electric field operator (local copy)
                g_ELocal = reducedFieldSILocal^2*bz.g_E;
                % evaluate power gained from the field
                power.field = gamma*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
                % evaluate correction due to the electric field term of the spatial growth model
                g_ELocal = reducedFieldSILocal*bz.g_fieldSpatialGrowth;
                correction = gamma*sum(eedfLocal.*(-g_ELocal(2:end)-g_ELocal(1:end-1)))*energyStep;
                power.field = power.field+correction;
    
                % diffusion contribution
                powerDiffusion = alphaRedEffLocal^2*gamma*energyStep/3*sum(energyCell(1:N).^2.*eedfLocal(1:N)./...
                    cellTotalCrossSection(1:N));
    
                % mobility contribution
                powerMobility = gamma*alphaRedEffLocal*(reducedFieldSILocal/6)*(...
                    energyCell(1)^2*eedfLocal(2)/cellTotalCrossSection(1) - ...
                    energyCell(N)^2*eedfLocal(N-1)/cellTotalCrossSection(N) + ...
                    sum(energyCell(2:N-1).^2.*(eedfLocal(3:N)-eedfLocal(1:N-2))./cellTotalCrossSection(2:N-1)));
    
                % power of spatial growth component
                power.eDensGrowth =  powerDiffusion + powerMobility;
        end
    else
        % elements of the electric field operator (local copy)
        g_ELocal = bz.workCond.reducedFieldSI^2*bz.g_E;
        % evaluate power gained from the field
        power.field = gamma*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
    end
    
    % power absorved per electron at unit gas density due to electron-electron collisions (to be revised)
    if bz.includeEECollisions
        A = bz.Afinal;
        B = bz.Bfinal;
        power.electronElectronNet = gamma*sum((A(1:N)-B(1:N)).*eedfLocal(1:N)')*energyStep^2;
        power.electronElectronGain = gamma*0.5*( sum((A(2:N-1)+B(3:N)-A(1:N-2)-B(2:N-1)).*eedfLocal(2:N-1)')+...
            (A(1)+B(2))*eedfLocal(1)-(A(N-1)+B(N))*eedfLocal(N) )*energyStep^2;
        power.electronElectronLoss = power.electronElectronNet - power.electronElectronGain;
    end
    
    % evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic collisions
    % loop over each gas in the mixture
    for gas = bz.gasArray
        gasName = gas.name;
        % initialize power balance information of this gas
        gasPower = struct('excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
            'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
            'ionizationIne', 0, 'attachmentIne', 0);
        % loop over each collision with the gas
        for collision = gas.collisionArray
            % collision type
            collType = collision.type;
            % avoid Effective or Elastic collisions and collisions which threshold is larger than the maximum energy
            if strcmp(collType, 'Effective') || strcmp(collType, 'Elastic') || collision.threshold > energyNode(end)
                continue;
            elseif strcmp(collision.type, 'Ionization') && ~strcmp(bz.ionCollOpType,'conservative')
                % evaluate cross section at cell positions
                cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
    
                threshold = collision.threshold;
                lmin=floor(threshold/energyStep);
    
                switch bz.ionCollOpType
                    case 'equalSharing'
                        ionizationIneAux = -gamma*collision.target.density*energyStep*(...
                            sum(energyCell(lmin:N).^2.*cellCrossSection(lmin:N).*eedfLocal(lmin:N)) + ...
                            2*energyCell(lmin+1)*sum(energyCell(2+lmin:2:N).*cellCrossSection(2+lmin:2:N).*...
                            eedfLocal(2+lmin:2:N))-2*sum(energyCell(2+lmin:2:N).^2.*cellCrossSection(2+lmin:2:N).*...
                            eedfLocal(2+lmin:2:N)));
    
                    case 'oneTakesAll'
                        ionizationIneAux = -gamma*collision.target.density*energyStep*energyCell(lmin)*...
                            sum(eedfLocal(lmin:N).*energyCell(lmin:N).*cellCrossSection(lmin:N));
    
                    case 'usingSDCS'
                        %calculation of the total (integrated) ionization cross section using the SDCS
                        TICS = zeros(size(eedfLocal));
                        W = gas.OPBParameter;
                        if isempty(W)
                            W = threshold;
                        end
                        auxArray = 1./(1+(energyCell(1:N)/W).^2);
                        for k=2:N
                            auxArray(k) = auxArray(k)+auxArray(k-1);
                            kmax = floor((k-lmin)/2);
                            if kmax>0
                                TICS(k) = TICS(k) + cellCrossSection(k)/(W*atan((energyCell(k)-threshold)/(2*W)))*auxArray(kmax);
                            end
                        end
                        TICS = energyStep*TICS;
                        ionizationIneAux = -gamma*collision.target.density*energyCell(lmin+1)*energyStep*...
                            sum(eedfLocal(1:N).*energyCell(1:N).*TICS(1:N));
                end
                gasPower.ionizationIne = gasPower.ionizationIne+ionizationIneAux;
                continue;
            elseif strcmp(collision.type, 'Attachment') && bz.includeNonConservativeAttachment
                % evaluate cross section at cell positions
                cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
                threshold = collision.threshold;
                lmin=floor(threshold/energyStep);
                gasPower.attachmentIne = gasPower.attachmentIne - gamma*collision.target.density*energyStep*...
                    sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).^2.*cellCrossSection(1+lmin:N));
                continue;
    
            end
            % switch to lower case because of aesthetical reasons
            collType = lower(collType);
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
            % evaluate departure cell
            lmin=floor(collision.threshold/energyStep);
            % add contribution to the power due to the inelastic collisions
            gasPower.([lower(collType) 'Ine']) = gasPower.([lower(collType) 'Ine']) - gamma*collision.target.density*...
                energyStep*energyNode(lmin+1)*sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
            % add contribution to the power due to the superelastic collisions
            if collision.isReverse
                statWeightRatio = collision.target.statisticalWeight/collision.productArray.statisticalWeight;
                gasPower.([lower(collType) 'Sup']) = gasPower.([lower(collType) 'Sup']) + gamma*statWeightRatio*...
                    collision.productArray.density*energyStep*energyNode(lmin+1)*sum(eedfLocal(1:N-lmin).*...
                    energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
            end
        end
        % evaluate net values (for each gas)
        gasPower.excitationNet = gasPower.excitationIne + gasPower.excitationSup;
        gasPower.vibrationalNet = gasPower.vibrationalIne + gasPower.vibrationalSup;
        gasPower.rotationalNet = gasPower.rotationalIne + gasPower.rotationalSup;
        gasPower.inelastic = gasPower.excitationIne + gasPower.vibrationalIne + gasPower.rotationalIne + ...
            gasPower.ionizationIne + gasPower.attachmentIne;
        gasPower.superelastic = gasPower.excitationSup + gasPower.vibrationalSup + gasPower.rotationalSup;
    
        % evaluate net values (for the gas mixture)
        power.excitationIne = power.excitationIne + gasPower.excitationIne;
        power.excitationSup = power.excitationSup + gasPower.excitationSup;
        power.vibrationalIne = power.vibrationalIne + gasPower.vibrationalIne;
        power.vibrationalSup = power.vibrationalSup + gasPower.vibrationalSup;
        power.rotationalIne = power.rotationalIne + gasPower.rotationalIne;
        power.rotationalSup = power.rotationalSup + gasPower.rotationalSup;
        power.ionizationIne = power.ionizationIne + gasPower.ionizationIne;
        power.attachmentIne = power.attachmentIne + gasPower.attachmentIne;
    
        % store power balance information of this gas in the main structure
        power.gases.(gasName) = gasPower;
    end
    power.excitationNet = power.excitationIne + power.excitationSup;
    power.vibrationalNet = power.vibrationalIne + power.vibrationalSup;
    power.rotationalNet = power.rotationalIne + power.rotationalSup;
    power.inelastic = power.excitationIne + power.vibrationalIne + power.rotationalIne + power.ionizationIne + ...
        power.attachmentIne;
    power.superelastic = power.excitationSup + power.vibrationalSup + power.rotationalSup;
    
    % evaluate power balance
    powerValues = [power.field power.elasticGain power.elasticLoss power.carGain power.carLoss ...
        power.excitationSup power.excitationIne power.vibrationalSup power.vibrationalIne ...
        power.rotationalSup power.rotationalIne power.eDensGrowth power.electronElectronGain ...
        power.electronElectronLoss];
    totalGain = 0;
    totalLoss = 0;
    for powerValue = powerValues
        if powerValue > 0
            totalGain = totalGain+powerValue;
        else
            totalLoss = totalLoss+powerValue;
        end
    end
    power.balance = power.field + power.elasticNet + power.carNet + power.inelastic + power.superelastic + ...
        power.eDensGrowth + power.electronElectronNet;
    power.relativeBalance = abs(power.balance)/totalGain;
    power.reference = totalGain;
    
    % store power balance information in the boltzmann properties
    bz.power = power;
    
    % check for errors in the power balance for final solution
    if checkPowerBalance && power.relativeBalance > bz.maxPowerBalanceRelError
        warning(sprintf(['Relative power balance greater than %e.\n' ...
            'Results may be wrong, please check input/output of the simulation'], bz.maxPowerBalanceRelError));
    end

end