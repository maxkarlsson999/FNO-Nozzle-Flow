function [data_solution,  branchType]=ExactSteadyStateSolution_Quasi1DNozzleFlow_Thesis(n_points, A, BC)
    global Cp Cv gamma R; 
    P0 = BC.P0_L ;      % Stagnation Pressure
    T0 = BC.T0_L ;      % Stagnation Temperature
    rho0 = P0/(R*T0);   % Stagnation Density
    P_e = BC.P_R  ;     % Exit pressure
    minA = min(A);      % Throat area
    A_e = A(end);       % Exit area
    iMin = find (A == minA);  
    iMin = iMin(1);    

    % Area-Mach number relation  with LHS = (A_e/ minA)^2
    area_mach_relation= @(M)    1/M^2 * (2/(gamma+1) * (1+ (gamma-1)/2 * M^2 ) )^((gamma+1)/(gamma-1)) ; 
    
    % Relation between Mach number before and after shock
    M_after_shock = @(M) sqrt((2 + (gamma - 1)*M^2)/(2*gamma*M^2 -(gamma - 1)));

    % Relation between pressure before and after shock
    P_after_shock = @(M, P) P + (2*P*gamma*(M^2 - 1))/(gamma + 1);

    % Ratio of unknown pressure and stagnation pressure 
    M_to_ratioP_P0 = @(M) ( 1 + (gamma-1)/2* M^2 )^( -gamma/(gamma-1) );
    
    % Mach Number from P/P0 ratio
    ratioP_P0_to_M = @(ratioP_P0)  sqrt( (ratioP_P0^(-(gamma-1)/gamma ) -1 )/(gamma-1)*2 );
    
    % Area-Mach number relation  with LHS = 0
    eqn__A_to_M = @(M) area_mach_relation(M) - ( A_e/ minA)^2 ;
    
    
    % Mach numbers for subsonic, shock at exit and supersonic cases
    M_e1 = fzero(eqn__A_to_M, [0.001, 0.9999]);
    M_e2 = M_after_shock(M_e1);
    M_e3 = fzero(eqn__A_to_M, [1, 100]);
    
    % Pressure for subsonic, shock at exit and supersonic cases
    P_e1 = M_to_ratioP_P0(M_e1)*P0
    P_e3 = M_to_ratioP_P0(M_e3)*P0
    P_e2 = P_after_shock(M_e3, P_e3)

    % Subsonic flow
    if  P_e >= P_e1
        branchType = 'subsonic';
        M_e = ratioP_P0_to_M(P_e/P0);
        minA_star = A_e/sqrt( area_mach_relation(M_e) );
        for i = 1: n_points
                A2 = ( A(i) / minA_star )^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
                M(i) = fzero(eqn__A_to_M, [0.0001, 0.9999]); 

                T(i)  = T0 ./ ( 1+ (gamma-1)/2 * M(i).^2 ); 
                P(i)  = P0./( (T0./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0 ./ ( (T0./T(i)).^(1/(gamma-1))   );             
        end
    end
    
    % Supersonic
    if P_e <= P_e2
       branchType = 'supersonic';
       for i = 1: n_points
            if( i<= iMin)
                A2 = ( A(i) / minA )^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
      
                M(i) = fzero( eqn__A_to_M, [0.0001, 1]); 
                T(i)  = T0 ./ ( 1+ (gamma-1)/2 * M(i).^2 ); 
                P(i)  = P0./( (T0./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0 ./ ( (T0./T(i)).^(1/(gamma-1))   );             
            else
                A2 = ( A(i) /  minA)^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
                M(i) = fzero( eqn__A_to_M , [1,100]); 

                T(i)  = T0 ./ ( 1+ (gamma-1)/2 * M(i).^2 ); 
                P(i)  = P0./( (T0./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0 ./ ( (T0./T(i)).^(1/(gamma-1))   ); 
            end
        end

    end
    
    % Shock in diverging sector
    if P_e < P_e1 && P_e > P_e2
        branchType = 'supersonicShock';
        A_1 = Calc_ShockPos_InNozzle(P_e, P0, A_e, minA);

        %
        PA_e__vs__P0minA_e  =  (P_e/P0) * ( A_e / minA ); 
        g_1 =  1/ (gamma-1); 
        M2_e  = - g_1   + sqrt( g_1^2 +  2*g_1 * ( 2/(gamma+1) ) ^ ( ( gamma+1)*g_1 ) * PA_e__vs__P0minA_e^(-2)    ) ; 
        M_e = sqrt( M2_e);
        P0_e = ( ( 1+ (gamma-1)/2* M_e^2 ) ^(gamma*g_1) ) * P_e; 
        
        %
        minA_e = P0*minA / P0_e; 
        T0_e = ( (P0_e * minA_e) / ( (P0*minA)/sqrt(T0) )  )^2 ; 
        rho0_e = P0_e / (R* T0_e ) ;
        %
        
        A_diverging =  A(iMin:end);
        iShock =  find ( A_diverging >= A_1 ); 
        iShock = iShock(1);

        %
        for i = 1: n_points
            if( i< iMin) % converge pipe
                
                A2 = ( A(i) / minA )^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
                M(i) = fzero( eqn__A_to_M, [0.0001, 1]); 
                    
                T(i)  = T0 ./ ( 1+ (gamma-1)/2 * M(i).^2 );
                P(i)  = P0./( (T0./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0 ./ ( (T0./T(i)).^(1/(gamma-1))   );             
            elseif( i< iShock+iMin)
                
                A2 = ( A(i) /  minA)^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
                M(i) = fzero( eqn__A_to_M , [1,100]); 

                T(i)  = T0 ./ ( 1+ (gamma-1)/2 * M(i).^2 ); 
                P(i)  = P0./( (T0./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0 ./ ( (T0./T(i)).^(1/(gamma-1))   ); 
            else
                A2 = ( A(i) / minA_e )^2;
                eqn__A_to_M = @(M) area_mach_relation(M)- A2;
                
                M(i) = fzero( eqn__A_to_M , [0.0001, 1]); 

                T(i)  = T0_e ./ ( 1+ (gamma-1)/2 * M(i).^2 ); 
                P(i)  = P0_e./( (T0_e./T(i)).^(gamma/(gamma-1))   );
                rho(i)= rho0_e ./ ( (T0_e./T(i)).^(1/(gamma-1))   ); 

            end
        end
    end




    data_solution.T   = T; 
    data_solution.P = P; 
    data_solution.rho = rho; 
    data_solution.u   = M .*  sqrt(gamma*R*T);
    data_solution.M   = M;
    data_solution.s  =  Cv*log( P./( (rho).^gamma ) );
    

end

% Determine shockPosition from an aribitary static pressure at the exit 
function [A_1]  = Calc_ShockPos_InNozzle(P_exit, P0, A_e, minA)
    global Cp Cv gamma R; 
    AnalyticalEqn__s2_s1__P02_P01  = @(M, P02_div_P01 ) Cp*log(  (1+2*gamma/(gamma+1) *(M^2-1) ) * ( 2+(gamma-1)*M^2)/ ( (gamma+1)*M^2 )   )  -R * log(1+ 2*gamma/(gamma+1) *(M^2-1) ) + R*log(P02_div_P01) ;  
    M_to_ratioA2                  = @(M)    1/M^2 * (2/(gamma+1) * (1+ (gamma-1)/2 * M^2 ) )^((gamma+1)/(gamma-1)) ; 
    
    % The following method are taken from page 216, chapter 5,
    %  from the book "Modern compressible flow, Third Edition", Johan, D. Anderson
    % The mass flow rate across nozzle
    %mdot = P0*minA/sqrt(T0) * sqrt( gamma/R* ( ( 2/(gamma+1) )^( (gamma+1)/(gamma-1) )  ) ); 
    PA_e__vs__P0minA_e  =  (P_exit/P0) * ( A_e / minA ); 
    g_1 =  1/ (gamma-1); 
    M2_e  = - g_1   + sqrt( g_1^2 +  2*g_1 * ( 2/(gamma+1) ) ^ ( ( gamma+1)*g_1 ) * PA_e__vs__P0minA_e^(-2)    ) ; 
    M_e = sqrt( M2_e);
    if( M_e <=1 )  
        P0_e = ( ( 1+ (gamma-1)/2* M_e^2 ) ^(gamma*g_1) ) * P_exit; 

        P02_div_P01 =  P0_e/P0;
        eqn_M__from__P02_P01  = @(M) AnalyticalEqn__s2_s1__P02_P01 ( M, P02_div_P01 );
        M_1 =  fzero(eqn_M__from__P02_P01 ,[1,100] ); 
        A_1 = sqrt( M_to_ratioA2(M_1) ) * minA; 
    else
        if( P_exit < 1000)
            A_1 = 999;
        else
            A_1 = -999;
        end
        
    end
    
end
