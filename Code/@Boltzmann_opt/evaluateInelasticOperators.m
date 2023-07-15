function evaluateInelasticOperators(bz)
    % evaluateInelasticOperators is in charge of the evaluation of the
    % elements of the inelastic and superelastic collision operator (it
    % does not include the contribution due to ionization nor attachment
    % collisions)
    %
    % For more information see the notes about the discretisation of the
    % boltzmann equation.
    
    % define local copies of variables used multiple
    eg = bz.energyGrid;
    energyNode = eg.node;
    energyCell = eg.cell;
    deps = eg.step;
    N = eg.cellNumber;
    ind_diag = 1:N+1:N*N;
    
    % allocate memory for the matrix
    M = zeros(N);
    
    % loop over each gas in the mixture
    for gas = bz.gasArray

        % loop over each collision with the gas
        for coll = gas.collisionArray
            
            % Local collision
            ctype = coll.type;
            threshold = coll.threshold;
            
            % avoid Effective, Elastic, Ionization and attachment
            % collisions, also avoid collisions which threshold is larger
            % than the maximum energy or smaller than the energy step
            if any(strcmp(ctype, {'Effective', 'Elastic', ...
                                  'Ionization', 'Attachment'}))
                continue;
            end
            if threshold > energyNode(end) || threshold < deps
                continue;
            end

            % evaluate numerical threshold
            nThresh = floor(threshold/deps);

            % evaluate cross section at cell positions
            nCS = coll.crossSection;
            cCS = 0.5*(nCS(1:end-1) + nCS(2:end));

            % store density of the target of the collision
            Tden = coll.target.density;

            % load matrix
    
            % evaluation of core inelastic/superelastic matrix elements
            eng_CS = energyCell.*cCS;

            % evaluation of inleastic matrix elements
            V = Tden * eng_CS;
            
            % fill "exits" in the inelastic matrix
            M(ind_diag) = M(ind_diag) - V;

            % fill "entrances" in the inelastic matrix
            ind = N*nThresh+1:N+1:N*N;
            M(ind) = M(ind) + V(1+nThresh:N);
            
            % Reversible
            if coll.isReverse
                
                % Ratio of degeneracies
                statWeightRatio = coll.target.statisticalWeight ./ ...
                                  coll.productArray.statisticalWeight;

                % Product Density
                Pden = coll.productArray.density;

                % evaluation of superelastic matrix elements
                V = statWeightRatio * Pden * eng_CS;

                % fill "exits" in the superelastic matrix
                ind = 1:N+1:N*(N-nThresh);
                M(ind) = M(ind) - V(1+nThresh:N);

                % fill "entrances" in the superelastic matrix
                ind = 1+nThresh:N+1:N*(N-nThresh);
                M(ind) = M(ind) + V(1+nThresh:N);
                
            end
            
        end

    end
    
    % save discrete matrix in the boltzmann object
    bz.inelasticMatrix = M;
    
end