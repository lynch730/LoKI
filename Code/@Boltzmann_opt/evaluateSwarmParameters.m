function swarmParam = evaluateSwarmParameters(bz)
    
    % save local copies of different constants and variables
    gamma = Constant.gamma;
    energyNode = bz.energyGrid.node;
    energyCell = bz.energyGrid.cell;
    energyStep = bz.energyGrid.step;
    eedfLocal = bz.eedf;
    WoN = bz.workCond.reducedExcFreqSI;
    
    % evaluate auxiliary total momentum transfer cross section function (see documentation)
    totalCrossSectionAux = bz.totalCrossSection;
    if strcmp(bz.eDensGrowthModel,'temporal') && ...
            ( bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment )
        totalCrossSectionAux(2:end) = totalCrossSectionAux(2:end) + (bz.CIEff/gamma)./sqrt(energyNode(2:end));
    end
    
    % initialize transport parameters structure
    swarmParam = struct('redDiffCoeff', [], 'redMobility', [], 'redMobilityHF', [], 'redDiffCoeffEnergy', [], ...
        'redMobilityEnergy', [], 'redTownsendCoeff', [], 'redAttCoeff', [], 'meanEnergy', [], 'characEnergy', [], ...
        'Te', [], 'driftVelocity', []);
    
    % evaluate reduced diffusion coefficient
    swarmParam.redDiffCoeff = (2*gamma/3)*energyStep*sum(energyCell.*eedfLocal./...
        (totalCrossSectionAux(1:end-1)+totalCrossSectionAux(2:end)));
    
    % evaluate reduced mobility (DC expression)
    swarmParam.redMobility = -gamma/3*sum(energyNode(2:end-1).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
        (totalCrossSectionAux(2:end-1)));
    
    % evaluate complex HF reduced mobility (in case the excitation frequency is not zero)
    if WoN ~= 0
        swarmParam.redMobilityHF = -gamma/3*sum(energyNode(2:end-1).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
            (totalCrossSectionAux(2:end-1)+(WoN/gamma)^2./(energyNode(2:end-1).*totalCrossSectionAux(2:end-1))));
        swarmParam.redMobilityHF = swarmParam.redMobilityHF + 1i*(1/3)*sum(sqrt(energyNode(2:end-1)).*...
            (WoN./totalCrossSectionAux(2:end-1)).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
            (totalCrossSectionAux(2:end-1)+(WoN/gamma)^2./(energyNode(2:end-1).*totalCrossSectionAux(2:end-1))));
    end
    
    % evaluate reduced energy diffusion coefficient
    swarmParam.redDiffCoeffEnergy = 2*gamma/3*energyStep*sum(energyCell.^2.*eedfLocal./...
        (totalCrossSectionAux(1:end-1)+totalCrossSectionAux(2:end)));
    
    % evaluate reduced energy mobility
    swarmParam.redMobilityEnergy = -gamma/3*sum(energyNode(2:end-1).^2.*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
        totalCrossSectionAux(2:end-1));
    
    % evaluate drift velocity
    if strcmp(bz.eDensGrowthModel,'spatial') && ...
            ( bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment )
        swarmParam.driftVelocity = - swarmParam.redDiffCoeff*bz.alphaRedEff + ...
            swarmParam.redMobility*bz.workCond.reducedFieldSI;
    else
        swarmParam.driftVelocity = swarmParam.redMobility*bz.workCond.reducedFieldSI;
    end
    
    % evaluate reduced Townsend coefficient
    totalIonRateCoeff = 0;
    for gas = bz.gasArray
        for collision = gas.collisionArray
            if strcmp(collision.type, 'Ionization')
                totalIonRateCoeff = totalIonRateCoeff + collision.target.density*collision.ineRateCoeff;
            end
        end
    end
    swarmParam.redTownsendCoeff = totalIonRateCoeff / swarmParam.driftVelocity;
    
    % evaluate reduced attachment coefficient
    totalAttRateCoeff = 0;
    for gas = bz.gasArray
        for collision = gas.collisionArray
            if strcmp(collision.type, 'Attachment')
                totalAttRateCoeff = totalAttRateCoeff + collision.target.density*collision.ineRateCoeff;
            end
        end
    end
    swarmParam.redAttCoeff = totalAttRateCoeff / swarmParam.driftVelocity;
    
    % evaluate mean energy
    swarmParam.meanEnergy = sum(energyCell(1:end).^(1.5).*eedfLocal(1:end))*energyStep;
    
    % evaluate characteristic energy
    swarmParam.characEnergy = swarmParam.redDiffCoeff/swarmParam.redMobility;
    
    % evaluate electron temperature
    swarmParam.Te = (2./3.)*swarmParam.meanEnergy;
    bz.workCond.update('electronTemperature', swarmParam.Te);
    
    % store swarm parameters information in the boltzmann properties
    bz.swarmParam = swarmParam;

end