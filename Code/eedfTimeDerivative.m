function dydt = eedfTimeDerivative(t, y, bz, cpv, ne_td_flag)
    % eedfTimeDerivative evaluates the time derivative of the EEDF at a
    % given time t, assuming that the EEDF at this precise time is given by
    % eedf.
    
    % Intermediates
    
    % separate variables into components
    eedf = y(1:end-1);    % electron energy distribution function
    ne = y(end);          % electron density (SI units)
    
    % local copy of fundamental constants
    persistent e;
    persistent e0;
    persistent gamma;
    if isempty(e)
        e = Constant.electronCharge;                % electron charge (SI units)
        gamma = Constant.gamma;                     % gamma parameter (sqrt(2e/me))
        e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
    end
    
    % local copy of other simulation constants and values
    persistent cellNumber;
    persistent energyStep;
    persistent energyCell;
    persistent energyNode;
    persistent N;
    persistent WoN;
    if isempty(cellNumber)
        cellNumber = bz.energyGrid.cellNumber;   % number of cells used in the energy grid
        energyStep = bz.energyGrid.step;         % energy step of the energy grid
        energyCell = bz.energyGrid.cell;         % energy at cells of the energy grid
        energyNode = bz.energyGrid.node;         % energy at nodes of the energy grid
        N = bz.workCond.gasDensity;              % gas density (SI units)
        WoN = bz.workCond.reducedExcFreqSI;      % reduced angular exitation frequency (SI units)
    end
    
    % evaluation of the reduced electric field (pulsed or constant value)
    if bz.isTimeDependent
        EoN = bz.pulseFunction(t,bz.pulseFunctionParameters)*1e-21;
    else
        EoN = bz.workCond.reducedFieldSI;
    end
    
    % renormalize EEDF
    eedf = eedf/sum(eedf.*sqrt(energyCell')*energyStep);
    
    % evaluate time idependent elements of the Boltzmann equation
    % (persistent for performance reasons)
    persistent matrix;
    persistent fieldMatrix;
    persistent nonConservativeMatrix;
    persistent CIEffIntegrand;
    persistent totalCrossSectionAux;
    persistent cellTotalCrossSection;
    persistent growthMatrixDiagElements;
    persistent growthMatrixSupElements;
    persistent growthMatrixInfElements;
    persistent g_extraFieldSpatialGrowth;
    persistent eeMatrixAuxA;
    persistent eeMatrixAuxB;
    if isempty(matrix)

        % evaluate basic boltzmann matrix (without field, ionization,
        % attachment or e-e collisions operators)
        matrix = bz.elasticMatrix + bz.CARMatrix + bz.inelasticMatrix;

        % save local copy of the regular field operator (without growth
        % models)
        fieldMatrix = bz.fieldMatrix;
    
        % include ionization operator (either conservative or not)
        nonConservativeMatrix = zeros(cellNumber);
        if bz.includeNonConservativeIonization
            matrix = matrix + bz.ionizationMatrix;
            nonConservativeMatrix = nonConservativeMatrix + bz.ionizationMatrix;
        else
            matrix = matrix + bz.ionizationConservativeMatrix;
        end
    
        % include attachment operator (either conservative or not)
        if bz.includeNonConservativeAttachment
            matrix = matrix + bz.attachmentMatrix;
            nonConservativeMatrix = nonConservativeMatrix + bz.attachmentMatrix;
        else
            matrix = matrix + bz.attachmentConservativeMatrix;
        end
    
        % evaluate time independent elements of the growth model operators
        % (in case they are activated)
        if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment

            % evaluate integrand to calculate the effective ionization rate
            CIEffIntegrand = sum(nonConservativeMatrix)';

            % local copies of total momentum transfer cross sections (at
            % nodes and cells)
            totalCrossSectionAux = bz.totalCrossSection;
            cellTotalCrossSection = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
            
            % choose growth model for the electron density (either spatial
            % or temporal)
            switch bz.eDensGrowthModel
                case 'temporal'

                    % evaluate matrix elements of the temporal growth
                    % operator (e-density time variation term, dn/dt)
                    growthMatrixDiagElements = -sqrt(energyCell)/gamma;

                case 'spatial'

                    % evaluate the diffusion and mobility components of the
                    % spatial growth operator
                    growthMatrixDiagElements = ...    % diffusion component
                        energyCell./(3*cellTotalCrossSection);

                    growthMatrixSupElements = ...     % mobility component (diag sup)
                        1/(6*energyStep)*[0 energyCell(1:cellNumber-1)./(cellTotalCrossSection(1:cellNumber-1))];
                    growthMatrixInfElements = ...     % mobility component (diag inf)
                        -1/(6*energyStep)*[energyCell(2:cellNumber)./(cellTotalCrossSection(2:cellNumber)) 0];
                    
                    % evaluate components of the extra electric field operator of the spatial growth model
                    g_extraFieldSpatialGrowth = energyNode/energyStep./(6*totalCrossSectionAux);
                    g_extraFieldSpatialGrowth(1) = 0;
                    g_extraFieldSpatialGrowth(end) = 0;

            end
        end
    
        % evaluate time idependent elements of the e-e collision operator
        % (in case that e-e collisions are activated)
        if bz.includeEECollisions
            
            % evaluating auxiliary matrix used in the evaluation of the ee
            % upflux vector without multiplicative constant
            eeMatrixAuxA = zeros(cellNumber);
            auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);

            % because of power conservation, terms on the last row
            % (cellNumber,:) and on the first column (:,1) are zero
            for k=1:cellNumber-1
                eeMatrixAuxA(k,2:k) = auxEnergyArray(2:k);
                eeMatrixAuxA(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
            end

            % detailed balance condition
            for k=1:cellNumber-1
                eeMatrixAuxA(k,2:cellNumber) = sqrt(eeMatrixAuxA(k,2:cellNumber).*eeMatrixAuxA(1:cellNumber-1,k+1)');
            end

            % evaluating auxiliary matrix used in the evaluation of the ee downflux vector without multiplicative constant
            eeMatrixAuxB = transpose(eeMatrixAuxA);

        end
    
    end
    
    % evaluate time derivative (time dependent) of each discrete component
    % of the eedf (except e-e collisions operator)
    if bz.includeNonConservativeIonization || bz.includeNonConservativeAttachment
        
        % calculation of the current effective ionization rate (needed for
        % both growth models, spatial and temporal)
        CIEff = gamma*energyStep*dot(eedf,CIEffIntegrand);
        
        % copy CIEff into boltzmann properties (needed for power balance
        % and swarm parameters calculation)
        bz.CIEff = CIEff;
        
        % choose growth model for the electron density (either spatial or
        % temporal)
        switch bz.eDensGrowthModel
            case 'temporal'

                % writing the modified total momentum transfer cross
                % section (total momentum transfer cross section plus the
                % ionization rate divided by the electron velocity) needed
                % for the reevaluation of the field operator
                totalCrossSectionMod = totalCrossSectionAux + (CIEff/gamma)./sqrt(energyNode);
                
                % writing the electric field operator matrix of the
                % temporal growth model (with totalCrossSectionMod)
                g_fieldTemporalGrowthAux = energyNode./( 3*energyStep^2*totalCrossSectionMod.* ...
                    (1+(WoN/gamma)^2./(energyNode.*totalCrossSectionMod.*totalCrossSectionMod)) );
                g_fieldTemporalGrowthAux(1) = 0;
                g_fieldTemporalGrowthAux(end) = 0;
                
                % copy field terms into boltzmann properties (needed for
                % power balance and swarm parameters calculation)
                bz.g_fieldTemporalGrowth = g_fieldTemporalGrowthAux*energyStep^2;
                
                % evaluate time derivative of each discrete component of
                % the eedf (case with temporal growth models)
                dfdt = N*gamma*( ( ... % multiplicative constant to obtain proper units of time
                    matrix*eedf + ...           % full basic boltzmann matrix
                    (EoN^2*( -g_fieldTemporalGrowthAux(1:cellNumber)-g_fieldTemporalGrowthAux(2:cellNumber+1)) + ... % time dependent diagonal component
                    CIEff*growthMatrixDiagElements)'.*eedf + ...
                    EoN^2*[ 0 g_fieldTemporalGrowthAux(2:cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...  % time dependent inf. diagonal component
                    EoN^2*[ g_fieldTemporalGrowthAux(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...      % time dependent sup. diagonal component
                    )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
            
            case 'spatial'

                % evaluation of gas density times diffusion coefficient and
                % mobility times the electric field
                ND = gamma*energyStep*dot(growthMatrixDiagElements,eedf);
                muE = -gamma*energyStep*EoN*dot(growthMatrixSupElements+growthMatrixInfElements,eedf);
                
                % calculation of the effective reduced first Townsend
                % coefficient
                if muE^2-4*CIEff*ND < 0
                    alphaRedEff = CIEff/muE;
                else
                    alphaRedEff = (muE - sqrt(muE^2-4*CIEff*ND))/(2*ND);
                end
                
                % copy alphaRedEff into boltzmann properties (needed for
                % power balance and swarm parameters calculation)
                bz.alphaRedEff = alphaRedEff;
                
                % copy field terms into boltzmann properties (needed for
                % power balance and swarm parameters calculation)
                bz.g_fieldSpatialGrowth = alphaRedEff*energyStep*g_extraFieldSpatialGrowth;
                
                % evaluate time derivative of each discrete component of
                % the eedf (case with spatial growth models)
                dfdt = N*gamma*( ( ...   % multiplicative constant to obtain proper units of time
                    matrix*eedf+...               % full basic boltzmann matrix (without field nor growth operators)
                    (EoN^2*fieldMatrix(1:cellNumber+1:cellNumber*cellNumber) - ...              % time dependent diagonal component
                    alphaRedEff*EoN*(g_extraFieldSpatialGrowth(1:cellNumber)-g_extraFieldSpatialGrowth(2:cellNumber+1)) + ...
                    alphaRedEff^2*growthMatrixDiagElements)'.*eedf + ...
                    [ 0 EoN^2*(fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) - ...          % time dependent inf. diagonal component
                    alphaRedEff*EoN*g_extraFieldSpatialGrowth(2:cellNumber) + ...
                    alphaRedEff*EoN*growthMatrixInfElements(1:cellNumber-1)) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
                    [ EoN^2*fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...  % time dependent sup. diagonal component
                    alphaRedEff*EoN*g_extraFieldSpatialGrowth(2:cellNumber) + ...
                    alphaRedEff*EoN*growthMatrixSupElements(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
                    )./sqrt(energyCell'));        % divide by the square root of the energy to obtain derivative of the EEDF
        end
    else
        % evaluate time derivative of each discrete component of the eedf
        % (case without growth models)
        dfdt = N*gamma*( ( ... % multiplicative constant to obtain proper units of time
            matrix*eedf + ...           % full basic boltzmann matrix (without field operator contribution)
            EoN^2*(fieldMatrix(1:cellNumber+1:cellNumber*cellNumber))'.*eedf + ...
            EoN^2*[ 0 fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
            EoN^2*[ fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
            )./sqrt(energyCell'));      % divide by the square root of the energy to obtain derivative of the EEDF
    end
    
    % evaluate e-e contribution to the time derivative (time dependent) of each discrete component of the eedf
    if bz.includeEECollisions
        
        % electron temperature in eV Te = (2/3)*meanEnergy
        Te = (2/3)*energyStep*dot(energyCell.^(3/2),eedf);
        
        % Coulomb logarithm
        logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));
        
        % multiplicative constant
        eeConstant = (ne/N)*(e^2/(8*pi*e0^2))*logC;
        
        % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
        A = (eeConstant/energyStep)*(eeMatrixAuxA*eedf);
        B = (eeConstant/energyStep)*(eeMatrixAuxB*eedf);
        
        % copy A and B into boltzmann properties (needed for power balance calculations)
        bz.Afinal = A;
        bz.Bfinal = B;
       
        % add contribution to time derivative of each discrete component of the eedf due to e-e collisions
        dfdt = dfdt + N*gamma*( ( ... % multiplicative constant to obtain proper units of time
            (-A(1:cellNumber)-B(1:cellNumber)).*eedf + ...                       % time dependent diagonal component
            [ 0; A(1:cellNumber-1) ].*[ 0; eedf(1:cellNumber-1) ] + ...          % time dependent inf. diagonal component
            [ B(2:cellNumber); 0 ].*[ eedf(2:cellNumber); 0 ] ...                % time dependent sup. diagonal component
            )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
    
    end
    
    % flush persistent memory for a new integration of the Boltzmann equation
    if cpv
        vars = whos;
        vars = vars([vars.persistent]);
        varName = {vars.name};
        clear(varName{:});
        dydt = [];
        return
    end
    
    % collect derivatives of the different variables
    if ne_td_flag
        dydt = [dfdt; ne*bz.workCond.gasDensity*CIEff];
    else
        dydt = [dfdt; 0];
    end

end
