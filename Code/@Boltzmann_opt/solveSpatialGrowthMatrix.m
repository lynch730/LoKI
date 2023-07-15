function eedf = solveSpatialGrowthMatrix(bz)
    
    % Local COpies
    eg = bz.energyGrid;
    energyCell = eg.cell;
    energyStep = eg.step;
    energyNode = eg.node;
    eedf = bz.eedf;
    N = eg.cellNumber;
    gamma = Constant.gamma;            % gamma parameter (sqrt(2e/me))
    EoN = bz.workCond.reducedFieldSI;  % reduced electric field (SI units)
    mixingParam = bz.mixingParameter;
    ind_d = 1:N+1:N*N;
    ind_dp = N+1:N+1:N*N;
    ind_dq = 2:N+1:N*N;
    
    % Initialize Matrices
    ionizationMatrixAux = bz.ionizationMatrix;
    attachmentMatrixAux = bz.attachmentMatrix;
    D = zeros(N);
    U = zeros(N);
    MFSG = zeros(N); % Spatial Growth
    sig_total = bz.totalCrossSection;
    cellTotalCrossSectionAux = 0.5*(sig_total(1:end-1) + sig_total(2:end));
    
    % writing of the Boltzmann matrix without the growth model and the
    % electron-electron collisional operator
    Mbase = bz.elasticMatrix + ...
            bz.CARMatrix + ...
            EoN^2*bz.fieldMatrix + ...
            bz.inelasticMatrix + ...
            bz.ionizationMatrix + ...
            bz.attachmentMatrix;
    
    % electron-electron collisions terms
    Mee = zeros(N);
    if bz.includeEECollisions
        alphaEE = bz.alphaee;
        auxA=bz.Aee;
        auxB=bz.Bee;
        A = (alphaEE/energyStep)*(auxA*eedf');
        B = (alphaEE/energyStep)*(auxB*eedf');
        Mee(1:N+1:N*N) = -(A(1:N)+B(1:N));
        Mee(ind_dp) = B(2:N);
        Mee(ind_dq) = A(1:N-1);
    end
    
    
    % evaluation of the effective ionization rate integrand
    integrandCI = gamma*energyStep*sum(ionizationMatrixAux+attachmentMatrixAux)';
    
    % evaluation of the diffusion and mobility components of the spatial
    % growth terms
    D0 = energyCell./(3*cellTotalCrossSectionAux);
    U0sup = EoN/(6*energyStep)*[0 energyCell(1:N-1)./(cellTotalCrossSectionAux(1:N-1))];
    U0inf = -EoN/(6*energyStep)*[energyCell(2:N)./(cellTotalCrossSectionAux(2:N)) 0];
    U0 = U0sup+U0inf;
    
    % evaluate effective ionization rate, reduced diffusion coefficient,
    % and mobility times the electric field
    CI = dot(eedf,integrandCI);
    ND = gamma*energyStep*dot(D0,eedf);
    muE = -gamma*energyStep*dot(U0,eedf);
    
    % evaluate reduced townsend coefficient (initial eedf guess may lead to
    % negative root argument, in which case, it should be calculated
    % assuming that there is no electron density gradient)
    if muE^2-4*CI*ND < 0
        alphaRedEffNew = CI/muE;
    else
        alphaRedEffNew = (muE - sqrt(muE^2-4*CI*ND))/(2*ND);
    end
    
    % initialize cycle counter
    iter = 0;
    max_iter = 400;
    while iter < max_iter 
    
        % writing of MatrixFieldSpatialGrowth which refers to the
        % additional electric field terms of the spatial growth model
        g_fieldSpatialGrowthAux = alphaRedEffNew*energyNode./(6*sig_total);
        g_fieldSpatialGrowthAux(1) = 0;
        g_fieldSpatialGrowthAux(end) = 0;
        for k=1:N
            MFSG(k,k) = -EoN*(g_fieldSpatialGrowthAux(k) - ...
                         g_fieldSpatialGrowthAux(k+1))/energyStep;
            if k>1
                MFSG(k,k-1) = -EoN*g_fieldSpatialGrowthAux(k)/energyStep;
            end
            if k<N
                MFSG(k,k+1) = EoN*g_fieldSpatialGrowthAux(k+1)/energyStep;
            end
        end
    
        % calculation of the diffusion and mobility matrices of the spatial
        % growth model
        D(ind_d) = alphaRedEffNew*alphaRedEffNew*D0;
        U(ind_dp) = alphaRedEffNew*U0sup(2:N);
        U(ind_dq) = alphaRedEffNew*U0inf(1:N-1);
    
        % writting of the Boltzmann matrix (expansion of the following
        % expression because of performance reasons) matrixAux =
        % 1.e20*(MatrixI + D + U + baseMatrix + Mee)
        M = Mbase;
        M(ind_d)  =  M(ind_d) + ( MFSG(ind_d)  + D(ind_d)  + Mee(ind_d)  );
        M(ind_dp) = M(ind_dp) + ( MFSG(ind_dp) + U(ind_dp) + Mee(ind_dp) );
        M(ind_dq) = M(ind_dq) + ( MFSG(ind_dq) + U(ind_dq) + Mee(ind_dq) );
        M = M * 1.0e20;

        % save previous solution
        eedfOld = eedf;
        alphaRedEffOld = alphaRedEffNew;
    
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(M, eg);
    
        % evaluate effective ionization rate, reduced diffusion
        % coefficient, and mobility times the electric field
        CI = dot(eedf,integrandCI);
        ND = gamma*energyStep*dot(D0,eedf);
        muE = -gamma*energyStep*dot(U0,eedf);
        
        % calculation of the new effective reduced first Townsend
        % coefficient
        if muE^2-4*CI*ND < 0
            alphaRedEffNew = CI/muE;
        else
            alphaRedEffNew = (muE - sqrt(muE^2 - 4*CI*ND))/(2*ND);
        end
    
        % evaluate convergence criteria
        if max(abs(eedf-eedfOld)./eedfOld) < bz.maxEedfRelError && ...
                (alphaRedEffNew==0 || (alphaRedEffNew-alphaRedEffOld)/alphaRedEffOld<1e-9)
            break;
        elseif iter == max_iter
            if ~bz.includeEECollisions
                warning('Spatial growth iterative scheme did not converge\n');
            end
            break;
        end
    
        % mixing of solutions
        alphaRedEffNew = alphaRedEffNew*mixingParam;
        alphaRedEffNew = alphaRedEffNew + (1-mixingParam)*alphaRedEffOld;
    
        % update cycle counter
        iter = iter + 1;
    
    end
    
    % copy of terms used in power balance and swarm parameters calculation
    bz.alphaRedEff = alphaRedEffOld;
    bz.g_fieldSpatialGrowth = g_fieldSpatialGrowthAux;
    bz.eedf = eedf;
    
    % copy of spatial growth model terms used on the electron-electron
    % collisions routine
    if bz.includeEECollisions
        bz.ionSpatialGrowthD = D;
        bz.ionSpatialGrowthU = U;
        bz.fieldMatrixSpatGrowth = MFSG;
    end
    
    % saving eedf
    bz.eedf = eedf;

end