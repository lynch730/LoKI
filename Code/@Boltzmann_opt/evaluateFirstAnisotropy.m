function evaluateFirstAnisotropy(bz)
    % evaluateFirstAnisotropy evaluates the value of the first anisotropy of the electron distribution function in
    % energy space in the framework of the two term approximation.
    
    % local copy of variables
    localEedf = bz.eedf;                             % electron energy distribution function (isotropic)
    energyCell = bz.energyGrid.cell;                 % values of energy at cell position (same as eedf)
    energyStep = bz.energyGrid.step;                 % energy step of the energy grid
    localTotalCrossSection = bz.totalCrossSection;   % total momentum transfer cross section
    gamma = Constant.gamma;                                 % gamma parameter (sqrt(2e/me))
    EoN = bz.workCond.reducedFieldSI;                % reduced electric field (SI units)
    WoN = bz.workCond.reducedExcFreqSI;              % reduced angular exitation frequency (SI units)
    
    % evaluate derivative of the eedf
    eedfDerivative = zeros(size(localEedf));
    eedfDerivative(1) = (localEedf(2)-localEedf(1))/energyStep;                     % 1st order forward approximation
    eedfDerivative(end) = (localEedf(end)-localEedf(end-1))/energyStep;             % 1st order backward approximation
    eedfDerivative(2:end-1) = (localEedf(3:end)-localEedf(1:end-2))/(energyStep*2); % 2nd order centered approximation
    
    % evaluate total momentum transfer cross section at cell positions
    totalCrossSectionCell = (localTotalCrossSection(1:end-1)+localTotalCrossSection(2:end))/2.0;
    
    % evaluate the first anisotropy
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        switch bz.eDensGrowthModel
            case 'temporal'
                totalCrossSectionCell = totalCrossSectionCell + (bz.CIEff/gamma)./(sqrt(energyCell));
                if WoN == 0
                    bz.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
                else
                    bz.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+(WoN/gamma)^2./...
                        (energyCell.*totalCrossSectionCell));
                end
            case 'spatial'
                bz.firstAnisotropy = -(bz.alphaRedEff.*localEedf+EoN.*eedfDerivative)./totalCrossSectionCell;
        end
    elseif WoN == 0
        bz.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
    else
        bz.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+(WoN/gamma)^2./...
            (energyCell.*totalCrossSectionCell));
    end

end