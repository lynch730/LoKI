function dydt = eedfTimeDerivative_opt(t, y, p)
    % eedfTimeDerivative evaluates the time derivative of the EEDF at a
    % given time t, assuming that the EEDF at this precise time is given by
    % eedf.
    
    % Current Field Strength
    EN = p.EoN(t);
    EN2 = EN * EN;
    
    % separate variables into components
    eedf = y(1:end-1);    % electron energy distribution function
    ne = y(end);          % electron density (SI units)
    
    % renormalize EEDF
    eedf = eedf / sum( eedf.*sqrt(p.eCell') * p.deps);
    
    % evaluate time derivative (time dependent) of each discrete component
    % of the eedf (except e-e collisions operator)
    if p.includeNonConservativeIonization || p.includeNonConservativeAttachment
        
        % calculation of the current effective ionization rate (needed for
        % both growth models, spatial and temporal)
        CIEff = p.gamma * p.deps * dot(eedf, p.CIEffIntegrand);
        
        % choose growth model for the electron density (either spatial or
        % temporal)
        switch p.eDensGrowthModel
            case 'temporal'
                
                % writing the modified total momentum transfer cross
                % section (total momentum transfer cross section plus the
                % ionization rate divided by the electron velocity) needed
                % for the reevaluation of the field operator
                tCSM = p.sig_total + (CIEff/p.gamma)./sqrt(p.eNode);
                
                % writing the electric field operator matrix of the
                % temporal growth model (with totalCrossSectionMod)
                g_fTGA = p.eNode./( 3*p.deps^2*tCSM.* ...
                    (1+(p.WoN/p.gamma)^2./(p.eNode.*tCSM.*tCSM)) );
                g_fTGA(1) = 0;
                g_fTGA(end) = 0;
                
                % evaluate w3q21`   time derivative of each discrete component of
                % the eedf (case with temporal growth models)
                %% MAIN COST FUNC
                
                g_fTGA = EN2 * g_fTGA;
                B = 1.0 ./ sqrt(p.eCell);
                
                M = p.M; 
                ind0 = 2:p.N+1:p.N*p.N;
                ind = 1:p.N+1:p.N*p.N;
                indp = p.N+1:p.N+1:p.N*p.N;
                M(ind0) = M(ind0) + g_fTGA(2:p.N) .* B(1:p.N-1);
                M(ind)  = M(ind)  + (- g_fTGA(1:p.N) - g_fTGA(2:p.N+1) + CIEff*p.M_Gdiag).*B;
                M(indp) = M(indp) + g_fTGA(2:p.N) .* B(2:p.N);

                M = M .* p.Nden .* p.gamma;
                dfdt = M * eedf;

            case 'spatial'

                % evaluation of gas density times diffusion coefficient and
                % mobility times the electric field
                ND = p.gamma*p.deps*dot(p.M_Gdiag,eedf);
                muE = -p.gamma*p.deps*EN*dot(p.M_GMSE+p.M_GMINFE,eedf);
                
                % calculation of the effective reduced first Townsend
                % coefficient
                if muE^2-4*CIEff*ND < 0
                    alphaRedEff = CIEff/muE;
                else
                    alphaRedEff = (muE - sqrt(muE^2-4*CIEff*ND))/(2*ND);
                end
                
                % evaluate time derivative of each discrete component of
                % the eedf (case with spatial growth models)
                dfdt = p.Nden*p.gamma*( ( ...   % multiplicative constant to obtain proper units of time
                    matrix*eedf+...               % full basic boltzmann matrix (without field nor growth operators)
                    (EN2*fieldMatrix(1:p.N+1:p.N*p.N) - ...              % time dependent diagonal component
                    alphaRedEff*EN*(g_extraFieldSpatialGrowth(1:p.N)-g_extraFieldSpatialGrowth(2:p.N+1)) + ...
                    alphaRedEff^2*p.M_Gdiag)'.*eedf + ...
                    [ 0 EN2*(fieldMatrix(2:p.N+1:p.N*p.N) - ...          % time dependent inf. diagonal component
                    alphaRedEff*EN*g_extraFieldSpatialGrowth(2:p.N) + ...
                    alphaRedEff*EN*p.M_GMINFE(1:p.N-1)) ]'.*[ 0; eedf(1:p.N-1) ] + ...
                    [ EN2*fieldMatrix(p.N+1:p.N+1:p.N*p.N) + ...  % time dependent sup. diagonal component
                    alphaRedEff*EN*g_extraFieldSpatialGrowth(2:p.N) + ...
                    alphaRedEff*EN*p.M_GMSE(2:p.N) 0 ]'.*[ eedf(2:p.N); 0 ] ...
                    )./sqrt(p.eCell'));        % divide by the square root of the energy to obtain derivative of the EEDF
        end
    else
        % evaluate time derivative of each discrete component of the eedf
        % (case without growth models)
        dfdt = p.Nden*p.gamma*( ( ... % multiplicative constant to obtain proper units of time
            matrix*eedf + ...           % full basic boltzmann matrix (without field operator contribution)
            EN2*(fieldMatrix(1:p.N+1:p.N*p.N))'.*eedf + ...
            EN2*[ 0 fieldMatrix(2:p.N+1:p.N*p.N) ]'.*[ 0; eedf(1:p.N-1) ] + ...
            EN2*[ fieldMatrix(p.N+1:p.N+1:p.N*p.N) 0 ]'.*[ eedf(2:p.N); 0 ] ...
            )./sqrt(p.eCell'));      % divide by the square root of the energy to obtain derivative of the EEDF
    end
    
    % evaluate e-e contribution to the time derivative (time dependent) of each discrete component of the eedf
    if p.includeEECollisions
        
        % electron temperature in eV Te = (2/3)*meanEnergy
        Te = (2/3)*p.deps*dot(p.eCell.^1.5,eedf);
        
        % Coulomb logarithm
        logC = log(12*pi*(p.eps0*Te/p.qe)^1.5/sqrt(ne));
        
        % multiplicative constant
        eeConstant = (ne/p.Nden)*(p.qe^2/(8*pi*p.eps0^2))*logC;
        
        % calculation of electron-electron collisions vectors of upflux (A)
        % and downflux (B)
        A = (eeConstant/p.deps)*(eeMatrixAuxA*eedf);
        B = (eeConstant/p.deps)*(eeMatrixAuxB*eedf);
        
        % add contribution to time derivative of each discrete component of
        % the eedf due to e-e collisions
        dfdt = dfdt + p.Nden*p.gamma*( ( ... % multiplicative constant to obtain proper units of time
            (-A(1:p.N)-B(1:p.N)).*eedf + ...                       % time dependent diagonal component
            [ 0; A(1:p.N-1) ].*[ 0; eedf(1:p.N-1) ] + ...          % time dependent inf. diagonal component
            [ B(2:p.N); 0 ].*[ eedf(2:p.N); 0 ] ...                % time dependent sup. diagonal component
            )./sqrt(p.eCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
        
    end
    
    % collect derivatives of the different variables
    if p.eDen_td_flag
        dydt = [dfdt; ne * p.ne_gasDensity * CIEff];
    else
        dydt = [dfdt; 0];
    end

end
