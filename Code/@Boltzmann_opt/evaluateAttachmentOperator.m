function evaluateAttachmentOperator(bz)
    % evaluateAttachmentOperator is in charge of the evaluation of the elements of the attachment
    % conservative or non-conservative collision operator
        
    % define local copies of variables used multiple times along the function
    eg = bz.energyGrid;
    energyNode = eg.node;
    energyCell = eg.cell;
    energyStep = eg.step;
    N = eg.cellNumber;
    Mcon = zeros(N);
    Mncon = zeros(N);
    thresh_flag = false;
    ind_diag = 1:N+1:N*N;

    % Loop Gases
    for gas = bz.gasArray

        % loop over each collision with the gas
        for coll = gas.collisionArray

            % Transition Energy Threshold
            thresh = coll.threshold;
            
            % Skip irrelevant
            if ~strcmp(coll.type, 'Attachment') || thresh > energyNode(end)
                continue
            end
            
            % set flag for non-conservative algorithms
            thresh_flag = true;
            
            % evaluate cross section at cell positions
            sig = 0.5*(coll.crossSection(1:end-1) + coll.crossSection(2:end));
            
            % local copy of target density
            Tden = coll.target.density;
            
            % evaluate numerical threshold
            NThresh = floor(thresh/energyStep);
            
            % evaluation of nonconservative attachment collisional operator
            % (always)
            Mncon(ind_diag) = Mncon(ind_diag) - Tden .* energyCell .* sig;
            
            % avoid writing of conservative operator if threshold is
            % smaller than the energy step
            if NThresh==0
                continue
            end

            % evaluation of conservative attachment collisional operator
            % (always needed)
            for k=1:N
                if (k<=N-NThresh)
                    Mcon(k,k+NThresh) = Mcon(k,k+NThresh) + ...
                                        Tden*energyCell(k+NThresh) * ...
                                        sig(k+NThresh);
                end
                Mcon(k,k) = Mcon(k,k) - Tden*energyCell(k)*sig(k);
            end

        end
    end
    
    % save matrixes in as properties of the boltzmann object
    bz.attachmentConservativeMatrix = Mcon;
    bz.attachmentMatrix = Mncon;
    bz.attachmentThresholdIsSmallerThanUmax = thresh_flag;
    
    if thresh_flag
        bz.includeNonConservativeAttachment = true;
    end

end