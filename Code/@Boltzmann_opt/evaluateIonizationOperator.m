function evaluateIonizationOperator(bz)
    % evaluateIonizationOperator is in charge of the evaluation of the elements of the ionization
    % conservative or non-conservative collision operator
    
    % define local copies of variables used multiple times along the function
    eg = bz.energyGrid;
    eNode = eg.node;
    eCell = eg.cell;
    deps = eg.step;
    N = eg.cellNumber;
    Mcon = zeros(N);
    Mncon = zeros(N);
    thresh_flag = false;
    ind_diag = 1:N+1:N*N;
    
    % Loop Gases
    for gas = bz.gasArray

        % loop over each collision with the gas
        for coll = gas.collisionArray

            % Transition energy thresdhold
            thresh = coll.threshold;
            
            % If Ionization
            if ~strcmp(coll.type, 'Ionization') || thresh > eNode(end)
                continue
            end
            
            % set flag for non-conservative algorithms
            thresh_flag = true;

            % evaluate cross section at cell positions
            sig = 0.5*(coll.crossSection(1:end-1) + coll.crossSection(2:end));
            
            % local copy of target density
            Tden = coll.target.density;

            % evaluate numerical threhold
            NThresh = floor(thresh/deps);

            % evaluate elements of the operator
            V = Tden.*eCell.*sig;

            % evaluation of nonconservative ionization collisional
            % operator (in case it is needed)
            switch bz.ionCollOpType
                case 'oneTakesAll'
                    
                    % fill "exits" in the ionNonConservativeMatrix
                    Mncon(ind_diag) = Mncon(ind_diag) - V;
                    
                    % fill "entrances" in the ionNonConservativeMatrix
                    % (scattered electrons)
                    ind = N*NThresh+1:N+1:N*N;
                    Mncon(ind) = Mncon(ind) + V(1+NThresh:N);

                    % fill "entrances" in the ionNonConservativeMatrix
                    % (secondary electrons)
                    Mncon(1,:) = Mncon(1,:) + V;

                case 'equalSharing'

                    % fill "exits" in the ionNonConservativeMatrix
                    Mncon(ind_diag) = Mncon(ind_diag) - V;

                    % fill "entrances" in the ionNonConservativeMatrix
                    % (both scattered and secondary electrons)
                    ind = N*(NThresh+1)+1:2*N+1:N*N;
                    indV = 2 + NThresh:2:2*floor((N-NThresh)/2)+NThresh;
                    Mncon(ind) = Mncon(ind) + 4*V(indV);
                    
                case 'usingSDCS'

                    W = gas.OPBParameter;
                    if isempty(W)
                        W = thresh;
                    end

                    % intermediate vars
                    aux1 = 1./(W.*atan((eCell-thresh)./(2*W)));
                    aux2 = 1./(1+(eCell./W).^2);
                    aux3 = deps * V .* aux1;
                    aux4 = cumsum(aux2);

                    % Loop Diagonal
                    for k=1:N

                        half = floor((k-NThresh)/2);
                        final = min([2*k+NThresh N]);
                        
                        % fill "exits" in the ionNonConservativeMatrix
                        if half>0
                            Mncon(k,k) = Mncon(k,k) - aux3(k)*aux4(half);
                        end

                        % fill "entrances" in the ionNonConservativeMatrix (both scattered and secondary electrons)
                        if k+NThresh+1<=N
                            ind = k+NThresh+1:final;
                            Mncon(k,ind) = Mncon(k,ind) + aux3(ind) .* ...
                                           aux2(1:final-k-NThresh);
                        end
                        
                        ind = (2*k+NThresh):N;
                        Mncon(k,ind) = Mncon(k,ind) + aux3(ind) .* aux2(k);
                        
                    end
            end

            % avoid writing of conservative collision operator if
            % threshold is smaller than the energy step
            if NThresh==0
                continue
            end

            % evaluation of conservative ionization collisional
            % operator (always needed) fill "exits" in the
            % ionConservativeMatrix
            Mcon(ind_diag) = Mcon(ind_diag) - V;
            
            % fill "entrances" in the ionConservativeMatrix
            ind = N*NThresh+1:N+1:N*N;
            Mcon(ind) = Mcon(ind) + V(1+NThresh:N);

        end
        
    end
    
    % save matrixes in as properties of the boltzmann object
    bz.ionizationConservativeMatrix = Mcon;
    bz.ionizationMatrix = Mncon;
    bz.ionizationThresholdIsSmallerThanUmax = thresh_flag;
    
    if ~strcmp(bz.ionCollOpType, 'conservative') && thresh_flag
        bz.includeNonConservativeIonization = true;
    end
    
end