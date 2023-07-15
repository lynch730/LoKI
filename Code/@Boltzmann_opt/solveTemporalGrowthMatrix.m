function eedf = solveTemporalGrowthMatrix(bz)
    
    gamma = Constant.gamma;              % gamma parameter (sqrt(2e/me))
    EoN = bz.workCond.reducedFieldSI;    % reduced electric field (SI units)
    WoN = bz.workCond.reducedExcFreqSI;  % reduced angular exitation frequency (SI units)
    
    eg = bz.energyGrid;
    eCell = eg.cell;
    deps = eg.step;
    N = eg.cellNumber;
    eNode = eg.node;
    eedf = bz.eedf;
    ind_d = 1:N+1:N*N;
    ind_p = N+1:N+1:N*N;
    ind_q = 2:N+1:N*N;
    mixingParam = bz.mixingParameter;
    
    % Init Matrices
    Mion = bz.ionizationMatrix;
    Matt = bz.attachmentMatrix;
    sig_total = bz.totalCrossSection;
    totalCSI = zeros(size(sig_total));
    g_TGA = zeros(size(bz.g_E));
    MFTG = zeros(N);
    Mgrowth = zeros(N);
    
    % writing of the Boltzmann matrix without the growth model and the
    % electron-electron collisional operator
    Mbase = bz.elasticMatrix + ...
                 bz.CARMatrix + ...
                 bz.inelasticMatrix + ...
                 Mion + Matt;
    
    % electron-electron collisions terms (not updated in the density growth
    % operator)
    Mee = zeros(N);
    if bz.includeEECollisions
        alphaEE = bz.alphaee;
        auxA=bz.Aee;
        auxB=bz.Bee;
        A = (alphaEE/deps)*(auxA*eedf');
        B = (alphaEE/deps)*(auxB*eedf');
        Mee(ind_d) = -(A(1:N)+B(1:N));
        Mee(ind_p) = B(2:N);
        Mee(ind_q) = A(1:N-1);
    end
    
    % evaluation of the effective ionization rate integrand
    integrandCI = gamma*deps*sum(Mion+Matt)';

    % evaluation of the initial value for the effective ionization rate
    CIEffNew = dot(eedf, integrandCI);
    
    % initialize cycle counter
    iter = 0;
    max_iter = 400;
    while iter < max_iter
        
        % writing the total cross section plus the ionization rate divided
        % by the electron velocity
        totalCSI(1)= sig_total(1);
        totalCSI(2:end) = sig_total(2:end) + (CIEffNew/gamma)./sqrt(eNode(2:end));
        
        % writing of the MatrixFieldTemporalGrowth which refers to the
        % electric field operator of the temporal growth model and
        % growthMatrix that refers to the time variation term (dn/dt) of
        % the temporal growth model
        g_TGA = eNode./(3*totalCSI.*(1+(WoN/gamma)^2./(eNode.*totalCSI.*totalCSI)));
        g_TGA(1) = 0;
        g_TGA(end) = 0;
        for k = 1:N
            MFTG(k,k) = -EoN^2*(g_TGA(k) + g_TGA(k+1))/deps^2;
            if k>1
                MFTG(k,k-1) = EoN^2*g_TGA(k)/deps^2;
            end
            if k<N
                MFTG(k,k+1) = EoN^2*g_TGA(k+1)/deps^2;
            end
            Mgrowth(k,k) = -(CIEffNew/gamma)*sqrt(eCell(k));
        end
    
        % writting of the Boltzmann matrix (expansion of the following
        % expression because of performance reasons)
        %         matrixAux =  1.e20*(growthMatrix + MatrixI + Mbase +
        %         Mee);
        M = Mbase;
        M(ind_d) = M(ind_d) + ( MFTG(ind_d) + Mgrowth(ind_d) + Mee(ind_d) );
        M(ind_p) = M(ind_p) + ( MFTG(ind_p) + Mee(ind_p) );
        M(ind_q) = M(ind_q) + ( MFTG(ind_q) + Mee(ind_q) );
        M = M * 1.0e20;

        % save previous solution
        eedfOld = eedf;
        CIEffOld = CIEffNew;
    
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(M, eg);
    
        % evaluate new effective ionization rate
        CIEffNew = dot(eedf,integrandCI);
    
        % evaluate convergence criteria
        if max(abs(eedf-eedfOld)./eedfOld)<bz.maxEedfRelError && ...
                (CIEffNew==0 || abs((CIEffNew-CIEffOld)/CIEffOld)<1e-9)
            break;
        elseif iter == max_iter
            if ~bz.includeEECollisions
                warning('Temporal growth iterative scheme did not converge');
            end
            break;
        end
    
        % mixing of solutions
        CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
    
        % update cycle counter
        iter = iter+1;
    end
    
    % copy of terms used in power balance and swarm parameters calculation
    bz.CIEff = CIEffOld;
    bz.g_fieldTemporalGrowth = g_TGA;
    bz.eedf = eedf;
    
    % copy of temporal growth model terms used on the electron-electron
    % collisions routine
    if bz.includeEECollisions
        bz.ionTemporalGrowth = Mgrowth;
        bz.fieldMatrixTempGrowth = MFTG;
    end
    
    % saving eedf
    bz.eedf = eedf;

end