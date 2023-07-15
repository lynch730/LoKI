function evaluateTotalAndElasticCrossSections(bz)
      
    % reset values to zero
    bz.totalCrossSection = zeros(1,bz.energyGrid.cellNumber+1);
    bz.elasticCrossSection = zeros(1,bz.energyGrid.cellNumber+1);

    % loop over each gas in the mixture
    for gas = bz.gasArray

        % avoid dummy gasses
        if isempty(gas.collisionArray)
            continue;
        end

        % evaluation of the mass ratio
        massRatio = Constant.electronMass / gas.mass;
        
        % loop over each collision with the gas
        for coll = gas.collisionArray

            coll_N = coll.target.density; % Temp copy
            
            % avoid effective collisions
            if strcmp(coll.type, 'Effective')
                continue;
            end
            
            % add collision cross section to the total momentum transfer
            %   cross section (also superelastic)
            sig = coll_N;
            if isempty(coll.momentumTransferCrossSection)
                sig = sig * coll.crossSection;
            else
                sig = sig * coll.momentumTransferCrossSection;
            end
            bz.totalCrossSection = bz.totalCrossSection + sig;
            
            % Do the same for reversible
            if coll.isReverse

                % Return Elastic and SuperElastic Cross Section
                [sE_CS, sE_MICS] = coll.superElasticCrossSection;

                % Density times cross section
                prod_N = coll.productArray.density;
                if isempty(sE_MICS)
                    prod_N = prod_N * sE_CS;
                else
                    prod_N = prod_N * sE_MICS;
                end
                bz.totalCrossSection = bz.totalCrossSection + prod_N;

            end

            % add elastic collision cross section to the total elastic
            %    cross section (weighted by the mass ratio)
            if strcmp(coll.type, 'Elastic')
                sig = massRatio*coll_N*coll.crossSection;
                bz.elasticCrossSection = bz.elasticCrossSection + sig;
                continue;
            end

        end
    end
    
end