
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This code demonstrates the implementaion of two numerical schemes (LAX 1st-order, Roe 1st-order) to solve 
%  the unsteady 2D Euler equations. 
%
%  The 2D solver has been setup to simulate the following problems for
%    (i)   VenturiNozzleFlow (included an exact quasi-1D analytical soluton for result comparison)
%
%
%  LTH course: MVKN70 (Advanced CFD ...)
%    Tutors: Rixin Yu
%            rixin.yu@energy.lth.se
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This code is for solve the unsteady 2D Euler equation:
%          d_t(U) + d_x(F) + d_y(G) = 0;
%   with
%         U = [rho,       rho*u,   rho*v ,  rho*E],  
%         F = [rho*u, rho*u^2+P, rho*u*v ,  rho*u*H],   
%         G = [rho*v, rho*u*v  , rho*v^2+P,  rho*u*H],   
%   with relations: 
%  (1) H = E + P/rho = Cp*T
%  (2) E = Cv*T + (u*u+v*v)/2
%  (3) P = rho*R*T 
%  where the above variables denote 
%   rho : density0
%   P   : pressure
%   u   : velocity
%   H   : total enthalpy (include kinetic energy)
%   E   : total energy   (include kinetic energy)
%   Cp  : specific heat capacity at constant pressure
%   Cv  : specific heat capacity at constant volume
%   R   : universal gas constant 
%
%   Note!!!! in this 2D implementaion, we choose a "core-set" "W"
%   formed by 4 unknowns in primitive form  
%          W =[rho,  u,  v ,  H],
%   any other quantities(such as P,T,.., etc.) can be expressed as 
%   a function of the three symbols inside W.  
%
%   After intial conditions being set, the solver enters a main loop performing explicit
%   time-advancement. The time-advacement is to update the values of three terms
%   inside the core-set W over all spatial mesh points from t^n to t^n+1, as
%   described in the following steps:
%    
%   Step 1:The core-set W are known at t^n (including t=0);
%   Step 2:  
%          Evalute all "relevant terms" according to certain numerical scheme using
%          the values of W given at t^n , such "relevant terms" includes (i) the 
%          four conservative variable in U, (ii) the four flux terms in F and (iii) etc.
%   Step 3:
%          Update U to  t^n+1 as: 
%
%           U(i, j; t^n+1)  =  SomeSchemeDependentFunction( W ;  i-1, i, i+1 ;  j-1, j, j+1 ; t^n) 
%   Step 4: 
%          Translate the four variables U(t^n+1) to W(t^n+1),  then enter the next loop
%   During the intermediate computaion, it is often convenient to switch between W , U , F ,  G.  
%


%   The nominal values used for air at 300 K are Cp = 1.00 kJ/kg.K, Cv = 0.718 kJ/kg.K,, and gamma = 1.4.
% The SI unit based on [kg, m, s , Kelvin] is used, note the unit of [joule=kg*(m/s)^2] [pascal=kg/m/s^2]

function [data, x, y, branch] = Euler2D_FunctionV3(P_L, P_R, T_0, r_nozzle)
    global Cp Cv gamma R; 
    global n_points;
    global x y x_1D;
    global jAll jF jo j1jo joj1 j_1jo joj_1 jojF jojF1 jFjo jF1jo;
    global CellVol  ProjectedFaceJacob ProjectedFaceDir;
    global Nosqueeze_CellVol Nosqueeze_ProjectedFaceJacob Nosqueeze_ProjectedFaceDir ;
    global dx; 
    global dt; 
    global data;
    global get_jAjF get_jFjA get_jF1jA get_jAjF1 get_jojo;
    global get_Jojo_half get_joJo_half  get_J_1jo_half get_joJ_1_half;
    global data_Anlytical_Quasi1D;
    global eps; % Roe's EntropyFix
    
    Cp    = 1000 ; % Do not Change
    gamma = 1.4;   % Do not Change
    Cv= Cp/gamma;  
    R = Cp - Cv; 
    %
    %%%%%
    
    
    BC.type = 'VenturiNozzleFlow'; 
    
    DomainSize_in_x =  5;  %[m]
    
    if contains( BC.type, 'VenturiNozzle' )
        
        % Three conditions needed to be prescribed!
        BC.P0_L = P_L;          % Total pressure    (inside a hypothetical zero-velocity gas reservoir)! [Pa]
        BC.T0_L = T_0;          % Total temperature (inside a hypothetical zero-velocity gas reservoir)! [K]
        BC.P_R  = P_R;          % Back/Exit pressure  [Pa]                              
        c = sqrt(  gamma*R* (BC.T0_L) );
        c_CharacteriticSoundSpeed =  c * 1.5; % Estimated max Mach number obtained inside Nozzle
    
        
        % nozzle geomerty&mesh
        n_points = [200 50]; % default total Number of grid cells 
        %n_points = [64 16]; % default total Number of grid cells 
        
        
          
    end
    
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Preperation for Mesh generation    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        
    x_1D      = linspace(-DomainSize_in_x/2,DomainSize_in_x/2 ,n_points(1));  x_1D=x_1D';
    x = zeros(n_points(1), n_points(2));   % 200x50
    y = zeros(n_points(1), n_points(2));   % 200x50
    for i = 1:n_points(1)
        x(i,:) = x_1D(i);  % copy the same x across all j
        y(i,:) = linspace(-r_nozzle(i), +r_nozzle(i), n_points(2));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Prepare mesh-related parameters for later caculation !!!! Begin  !!!!
    %              (read lecture 6, slides 5, 6 and 7)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %  1: 2D mesh "indexes"
    %
    
    jAll= { 1:n_points(1), 1:n_points(2) }; % the j-index used for all cells
    jF= { 1:n_points(1)-1, 1:n_points(2)-1 }; % the j-index used for F_jhalf 
    jo ={ 2:n_points(1)-1, 2:n_points(2)-1}; % The j-indexs used to refer the "interior" mesh points
    j1jo = {jo{1}+1 , jo{2} } ; j_1jo = {jo{1}-1 , jo{2} };
    joj1 = {jo{1} , jo{2}+1 } ; joj_1 = {jo{1} , jo{2}-1 };
    jojF = { jo{1}, jF{2} }   ; jojF1 = {jo{1} , jF{2}+1}; 
    jFjo = { jF{1}, jo{2} }   ; jF1jo = { jF{1}+1, jo{2}  }; 
    
    jFjA = {jF{1}, jAll{2} } ; jAjF = { jAll{1} , jF{2} };
    jF1jA = {jF{1}+1, jAll{2} } ; jAjF1 = { jAll{1} , jF{2}+1 } ;
    
    %
    %  2: Edge-mesh as the "dual" of cell-center-mesh
    %
    
    Edges.x( 1+jF{1}, 1+jF{2}  ) = ( x( jF{1},jF{2} ) + x( jF{1},jF{2}+1 ) +  x( jF{1}+1,jF{2} )  + x( jF{1}+1,jF{2}+1) )/4  ;
    Edges.y( 1+jF{1}, 1+jF{2}  ) = ( y( jF{1},jF{2} ) + y( jF{1},jF{2}+1 ) +  y( jF{1}+1,jF{2} )  + y( jF{1}+1,jF{2}+1) )/4  ;
    Edges.x( 1, : ) = 2*Edges.x( 2,:)-Edges.x( 3,:); Edges.x(end+1,:)= 2*Edges.x(end,:)-Edges.x(end-1,:);
    Edges.y( 1, : ) = 2*Edges.y( 2,:)-Edges.y( 3,:); Edges.y(end+1,:)= 2*Edges.y(end,:)-Edges.y(end-1,:); 
    Edges.x(:,    1)= 2*Edges.x(:, 2)- Edges.x(:,    3);Edges.x( :,end+1)=2*Edges.x(:,end)- Edges.x(:,end-1);
    Edges.y(:,    1)= 2*Edges.y(:,2 )-Edges.y(:,    3); Edges.y(:,end+1)= 2*Edges.y(:,end)-Edges.y(:,end-1)  ; 
    
    % Uncomment this to view your mesh
    %surf(Edges.x, Edges.y,Edges.y*0-1); daspect([1 1 1]); view(2); hold on; plot(Edges.x,Edges.y, 'xb'); hold on; plot(x,y, 'or'); hold on; return;
    
    %
    % 3:  Volume of mesh cells caculated using matlab function "polyarea" 
    %
    for jx = 1: n_points(1)
    for jy = 1: n_points(2)
     CellVol(jx,jy) = polyarea ( [ Edges.x(jx,jy), Edges.x(jx,jy+1), Edges.x(jx+1,jy+1), Edges.x(jx+1,jy)],  [ Edges.y(jx,jy), Edges.y(jx,jy+1), Edges.y(jx+1,jy+1), Edges.y(jx+1,jy)] );    
    end
    end
    % Uncomment this to plot the cell volume
    %mesh(x, y , CellVol) ;daspect([1 1 1]); view(2);pause(.2);return;
    
    Edges.jacob11 = Edges.x ( 2:end, 2:end)-Edges.x(1:end-1, 2:end )  ; 
    Edges.jacob12 = Edges.y ( 2:end, 2:end)-Edges.y(1:end-1, 2:end )  ; 
    Edges.jacob21 = Edges.x ( 2:end, 2:end)-Edges.x(2:end, 1:end-1)  ; 
    Edges.jacob22 = Edges.y ( 2:end, 2:end)-Edges.y(2:end, 1:end-1) ; 
    %
    % 4 : projected face jocobain (read slide 7 in lecture 6)
    %
    ProjectedFaceJacob{1,1} =   Edges.jacob22 ; 
    ProjectedFaceJacob{1,2} =  -Edges.jacob21 ; 
    ProjectedFaceJacob{2,1} =  -Edges.jacob12 ;
    ProjectedFaceJacob{2,2} =   Edges.jacob11 ;
    
    k = 1;
    for k2 = 1:2
       ProjectedFaceDir{k,k2} = sign( ProjectedFaceJacob{k,k2}(1:end-1,:)   );  % for determining the "upwinding" direction in a 2D cases (required for Roe scheme)
    end
    k = 2;
    for k2 = 1:2
       ProjectedFaceDir{k,k2} = sign( ProjectedFaceJacob{k,k2}(:,1:end-1)   );  % for determining the "upwinding" direction in a 2D cases (required for Roe scheme)
    end
    %
    % 5: Coordinate difference between neghbouring cell-centers (used for lax scheme & boundary conditions)
    % 
    dx{1,1} = permute(  repmat( x(jF1jA{:})-x( jFjA{:} )  , [1 1 4]) , [3 1 2] ); 
    dx{1,2} = permute(  repmat( y(jF1jA{:})-y( jFjA{:} )  , [1 1 4]) , [3 1 2] );
    dx{2,1} = permute(  repmat( x(jAjF1{:})-x( jAjF{:} )  , [1 1 4]) , [3 1 2] );
    dx{2,2} = permute(  repmat( y(jAjF1{:})-y( jAjF{:} )  , [1 1 4]) , [3 1 2] );
    %
    % For Computational speed-up (perfomance enhancing tricks only for matlab)
    %
    
    Nosqueeze_CellVol                      = permute(  repmat( CellVol  , [1 1 4]) , [3 1 2] ); 
    for k = 1:2
    for k2 = 1:2
       Nosqueeze_ProjectedFaceJacob{k,k2} = permute(  repmat( ProjectedFaceJacob{k,k2}  , [1 1 4]) , [3 1 2] ); 
       Nosqueeze_ProjectedFaceDir{k,k2}   = permute(  repmat( ProjectedFaceDir{k,k2}  , [1 1 4]) , [3 1 2] )  ;
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Prepare mesh-related parameters for later caculation: !!!! Ready  !!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TotalSimulationTime = DomainSize_in_x/c_CharacteriticSoundSpeed *100 ; %Change the total time of the simulation here if you want longer time of simulations

    CFL  = 0.3;
    
    delta_x =    sqrt(Edges.jacob11.^2 +Edges.jacob12.^2 ) ;
    delta_y =    sqrt(Edges.jacob21.^2+ Edges.jacob22.^2 ) ; % minium cell distance
    
    min_dx = min( min( min( sqrt( Edges.jacob11.^2 +Edges.jacob12.^2 ) ) ) ,  min( min( sqrt(Edges.jacob21.^2+ Edges.jacob22.^2 ) ) )); % minimum cell distance
     
    dt = CFL * min_dx /c_CharacteriticSoundSpeed;  % time-step size
    
    Num_T = round( TotalSimulationTime/dt ); % total number of time stepping
   
    Prepare_and_AllocateMemoryForComputionalSpeedUp(); % This is critical for performance
    %
    
    if contains( BC.type, 'VenturiNozzle'  ) 
        
        % This case has a quasi-1D anlytical steady-state solution, here we
        % obtain the analytical soluton to be compared with the later 2D results.
        area = pi*(r_nozzle).^2;
        [data_Anlytical_Quasi1D, branchType] =ExactSteadyStateSolution_Quasi1DNozzleFlow_Max(n_points(1), r_nozzle , BC); 
        branch = branchType;
        %plotNozzle(x_1D,r_nozzle); % show the analytical prediction from qusai-1D analysis 
        %pause(3);  % after 3 seconds see the results, we proceed! 
        %close all;
        
        %for j =1 : n_points(2)
        %    data.rho(:,j) = data_Anlytical_Quasi1D.rho(:);        data.u(:,j) = data_Anlytical_Quasi1D.u(:);        data.T(:,j) = data_Anlytical_Quasi1D.T(:);
        %end
        %data_Anlytical_Quasi1D.CurrentTime = '100';
             
    end
    
    %AllResults = {};  
    
    %%%%%%%%%%%%%%%%%%%%%%
    NameOfAvaiableSchemes={                          ...
         'LAX Scheme:1st order'                     ...
         'Roe Scheme:1st order EntropyFix',        ...
         'Roe Scheme:1st order',        ...
         } ;   
    tic
    
    i_Selected_Scheme = 2; 
    %for i_Selected_Scheme = 2
        str_NameOfNumericalScheme =  NameOfAvaiableSchemes{i_Selected_Scheme}; 
        
        if contains( str_NameOfNumericalScheme, 'EntropyFix')
            eps = 100;   
        else
            eps = 0;   
        end
       
        data = SetInitialConditions(data,BC); 
        
        % plot prepartion
        fig = figure(i_Selected_Scheme);   set(fig, 'DefaultLineLineWidth', 0.5);   set(fig, 'units','normalized','outerposition',[0 0 1 1]);
        data.CurrentTime = 0;
        tol = 1e-10;
        converged = false;
        for TimeStep=1: Num_T
             U_old = data.U;
    
            if contains( str_NameOfNumericalScheme, 'LAX Scheme:1st order')    
               % A reminder:  jF = 1: n_points-1; is the j-index used for F_jhalf 
              data.AreaFluxA_jhalf(:,:,:) = ...
                    (get_jFjA(data.F) + get_jF1jA(data.F) - ( get_jF1jA(data.U) - get_jFjA(data.U) ).* dx{1,1}/dt/2 )/2 .* get_jFjA( Nosqueeze_ProjectedFaceJacob{1,1} ) + ...
                    (get_jFjA(data.G) + get_jF1jA(data.G) - ( get_jF1jA(data.U) - get_jFjA(data.U) ).* dx{1,2}/dt/2 )/2 .* get_jFjA( Nosqueeze_ProjectedFaceJacob{1,2} ) ;
              data.AreaFluxB_jhalf(:,:,:) = ...
                    (get_jAjF(data.F) + get_jAjF1(data.F) - ( get_jAjF1(data.U) - get_jAjF(data.U) ).* dx{2,1} /dt/2 )/2.* get_jAjF( Nosqueeze_ProjectedFaceJacob{2,1} ) + ... 
                    (get_jAjF(data.G) + get_jAjF1(data.G) - ( get_jAjF1(data.U) - get_jAjF(data.U) ).* dx{2,2} /dt/2 )/2.* get_jAjF( Nosqueeze_ProjectedFaceJacob{2,2} ) ;
    
    
            elseif contains( str_NameOfNumericalScheme, 'Roe' )  % Schemes based on Roe's linearized-Riemann-solver 
    
               [data.AreaFluxA_jhalf, data.AreaFluxB_jhalf ]= Update_Roe_Flux(data, 'piecewise-constant' ); 
    
            else
               errordlg( [  'not implemented ' str_NameOfNumericalScheme] ); 
               return;
            end 
           
            
            % The general "finite-volume" form 
            data.U(:,2:end-1,2:end-1)= get_jojo(data.U) ...
                 - dt* (  get_Jojo_half( data.AreaFluxA_jhalf )-get_J_1jo_half(data.AreaFluxA_jhalf)  ...
                          + get_joJo_half( data.AreaFluxB_jhalf )- get_joJ_1_half(data.AreaFluxB_jhalf)  ...
                        ) ./ Nosqueeze_CellVol(:, 2:end-1,2:end-1 ); 
            data = U_to_W(jo,data);
            data=Set_BoundaryCondition(data,BC,dt,dx);
            data=Update_PTEsM_AfterReset_W(data);           
            data.CurrentTime = data.CurrentTime + dt;
            dU = data.U - U_old;
            residual = max(abs(dU(:))) / max(abs(U_old(:)));

            if residual < tol
                disp(['Converged after ', num2str(TimeStep), ' iterations, diff=', num2str(residual)]);
                converged = true;
                break;
            end
            %------
            if mod(TimeStep , 1 ) == 20 
               c = sqrt( gamma*R)*sqrt(data.T)  ;
               dt = 0.6 * min ( min( min( delta_x ./ (c + abs(data.u) ) )),  min(min( delta_y ./ (c + abs(data.v)) ) )  ) ;
            end
            %------
            
            
            %  Plot every 100 time step togher with analytical shock-tube solutions 
            if mod(TimeStep , 100 ) == 1 
                
               ElapsedTime=toc; strElapsedTime = sprintf('(Elapsed computer time: %ds)',round(ElapsedTime) );
               
               pause_plot_2D(x, y, data, [str_NameOfNumericalScheme strElapsedTime  ] ,BC );
              
               pause(.001)
               fprintf('Step %d, diff = %.3e\n', TimeStep, residual);
            end
            
        end % end loop of time-advancement 
       
        
        %AllResults{end+1} = { str_NameOfNumericalScheme, data}; 
         
    %end  % end loop of all two aviable schemes
   
    return;
   
end

% change here if you want other boundary conditions
function [data]= Set_BoundaryCondition(data, BC,dt,dx)
    if contains( BC.type, 'Nozzle')
        data=BoundaryCondition_Nozzle_LeftSubsonicInlet_RightOutlet(data, BC);
        data=BoundaryCondition_Wall(data,dt,1,1); % bottom wall;
        data=BoundaryCondition_Wall(data,dt,1,2); % top wall;

    else
        % all zero-gradient
        data=BoundaryCondition_ZeroGradient(data,1,1); 
        data=BoundaryCondition_ZeroGradient(data,1,2); 
        data=BoundaryCondition_ZeroGradient(data,2,1); 
        data=BoundaryCondition_ZeroGradient(data,2,2); 
    end
end

function [data] = Update_PTEsM_AfterReset_W(data)
    global Cp Cv gamma R; 
    data.T  = ( (data.H) - 0.5* ( (data.u).^2+(data.v).^2)  )/Cp ; 
    data.E  = Cv* (data.T) + 0.5* ( (data.u).^2+(data.v).^2 )  ;       % total energy
    data.P  = R* (data.rho) .* (data.T); 
    data.s  = Cv*log( data.P  ./  (data.rho).^gamma    ); 
    
    %data.M  = ( data.u ) ./ sqrt( gamma * R * data.T );  
   
    data.U(1,:,:) =  data.rho;                
    data.U(2,:,:) =  data.rho.*data.u; 
    data.U(3,:,:) =  data.rho.*data.v; 
    data.U(4,:,:) =  data.rho.*data.E; 
    
    rhou = data.rho .* data.u;
    data.F(1,:,:)=rhou; 
    data.F(2,:,:)=rhou .* data.u+data.P ;  
    data.F(3,:,:)=rhou .* data.v;  
    data.F(4,:,:)=rhou .* data.H;  
    
    rhov = data.rho .* data.v;
    data.G(1,:,:)=  rhov                        ;   
    data.G(2,:,:)=  rhov .* data.u;   
    data.G(3,:,:)=  rhov .* data.v+data.P ;
    data.G(4,:,:)=  rhov .* data.H            ;

end

% Translate the inputting U to the 3 primitive variable [rho,u,H] inside
% core-set W, then write to the correponding spatial memories carried in 
% the structure "data" at the location of j cell 
function data= U_to_W(j,data)
  global Cp Cv gamma R; 
  global jAll jF jo j1jo joj1 j_1jo joj_1 jojF jojF1 jFjo jF1jo;
  if  isequal( j, jo)
     data.rho(2:end-1, 2:end-1) = squeeze( data.U(1, 2:end-1, 2:end-1) ); 
     data.u( 2:end-1, 2:end-1 )   = squeeze( data.U(2, 2:end-1, 2:end-1 ) ./ data.U(1,  2:end-1, 2:end-1 ) ); 
     data.v( 2:end-1, 2:end-1 )   = squeeze( data.U(3, 2:end-1, 2:end-1) ./ data.U(1,  2:end-1, 2:end-1) ); 
     data.H( 2:end-1, 2:end-1 )   = squeeze( gamma*data.U(4,   2:end-1, 2:end-1 ) ./ data.U(1, 2:end-1, 2:end-1 )) - (gamma-1)*0.5* ( ( data.u(2:end-1, 2:end-1)).^2 + (data.v( 2:end-1, 2:end-1) ) .^2 ) ; 
  end  
   
end

function data_LR= U_to_data_LR(data_LR)
    global Cp Cv gamma R;
    data_LR.rho = squeeze( data_LR.U(1,:,:) ); 
    data_LR.u   = squeeze( data_LR.U(2,:,:) ./ data_LR.U(1,:,:) ) ; 
    data_LR.v   = squeeze( data_LR.U(3,:,:) ./ data_LR.U(1,:,:) ) ; 
    data_LR.H   = gamma*squeeze(data_LR.U(4,:,:)./ data_LR.U(1,:,:))  - (gamma-1)*0.5* (  (data_LR.u ).^2 +  (data_LR.v).^2 ) ; 


    data_LR.T  = ( (data_LR.H) - 0.5* (  (data_LR.u).^2  + (data_LR.v).^2 )  )/Cp ; 
    data_LR.P  = R* (data_LR.rho) .* (data_LR.T); 

    rhou = data_LR.rho .* data_LR.u;
    data_LR.F(1,:,:)=rhou; 
    data_LR.F(2,:,:)=rhou .* data_LR.u+data_LR.P ;  
    data_LR.F(3,:,:)=rhou .* data_LR.v;  
    data_LR.F(4,:,:)=rhou .* data_LR.H;  

    rhov = data_LR.rho .* data_LR.v;
    data_LR.G(1,:,:)=  rhov                        ;   
    data_LR.G(2,:,:)=  rhov .* data_LR.u;   
    data_LR.G(3,:,:)=  rhov .* data_LR.v+data_LR.P ;
    data_LR.G(4,:,:)=  rhov .* data_LR.H            ;
end

function Roe = get_invYdV__given_data_LR(Roe,data_L,data_R )
   dV{1} =  data_R.rho-data_L.rho   ; 
   dV{2} =  data_R.u  -data_L.u     ;
   dV{3} =  data_R.v  -data_L.v     ;
   dV{4} =  data_R.P  -data_L.P     ;

   % for eigendecompose(F) =   abs(lambda) .*     alpha (i.e. [pi]*inv[Y]*d_V)  
   k =1;
   Roe.invYdV_F{1} =                       -Roe.Ave.c.*Roe.Ave.rho.* dV{2}                      + dV{4} ; 
   Roe.invYdV_F{2} =  Roe.Ave.c.^2.*dV{1}                                                       - dV{4} ; 
   Roe.invYdV_F{3} =                                                      2*Roe.Ave.c.^2.*dV{3}         ;  
   Roe.invYdV_F{4} =                        Roe.Ave.c.*Roe.Ave.rho.* dV{2}                      + dV{4} ;

   % for  decompose(G), based on decompose(F),  just  switch u and v 
   k = 2;
   Roe.invYdV_G{1} =                       -Roe.Ave.c.*Roe.Ave.rho.* dV{3}                      + dV{4} ; 
   Roe.invYdV_G{3} =  Roe.Ave.c.^2.*dV{1}                                                       - dV{4} ; 
   Roe.invYdV_G{2} =                                                      2*Roe.Ave.c.^2.*dV{2}         ;  
   Roe.invYdV_G{4} =                        Roe.Ave.c.*Roe.Ave.rho.* dV{3}                      + dV{4} ;
end

function abs_uc_fixed = EntropyFix( abs_uc )
   global eps;
   abs_uc_fixed = abs_uc;    
   % Roe entropy fix (it will be disabled if eps is set too small)
   % eps=100 ; %to enable entropy fix, set a large value for eps, say 100 
   if( eps > 0.1) 
      abs_uc_fixed(abs_uc < eps) = ( abs_uc(abs_uc < eps).^2 + eps^2 ) / eps/2 ;         
   end
end

% [RoeAve_X]* ( abs[RoeAve_Lambda]* [PI] )
function [Roe]  = LeftMultiplyBy_RoeAve_X_AbsLambda_PI (Roe,  k  )
   global ProjectedFaceDir;
   c22 = 2*Roe.Ave.c.^2;
   k2 = 1;
   A = Roe.invYdV_F;
   abs__u_c  = EntropyFix( abs(Roe.Ave.u-Roe.Ave.c ) ).* ProjectedFaceDir{k,k2} ;
   abs__u    = abs(Roe.Ave.u)                        .* ProjectedFaceDir{k,k2}  ; 
   abs__uc   = EntropyFix(  abs(Roe.Ave.u+Roe.Ave.c) ).* ProjectedFaceDir{k,k2}   ; 
   Roe.XLAlpha_F(1,:,:) =  (abs__u_c                               .*A{1} + 2*abs__u                           .*A{2}+                                       abs__uc                              .*A{4})./c22; 
   Roe.XLAlpha_F(2,:,:) =  (abs__u_c.*(Roe.Ave.u-Roe.Ave.c )       .*A{1} + 2*abs__u.*Roe.Ave.u                 .*A{2}+                                       abs__uc.*(Roe.Ave.u+Roe.Ave.c)         .*A{4})./c22; 
   Roe.XLAlpha_F(3,:,:) =  (abs__u_c.*Roe.Ave.v                    .*A{1} + 2*abs__u.*Roe.Ave.v                 .*A{2}+ abs__u.*Roe.Ave.rho          .*A{3}+   abs__uc.*Roe.Ave.v                    .*A{4})./c22;
   Roe.XLAlpha_F(4,:,:) =  (abs__u_c.*(Roe.Ave.H-Roe.Ave.u.*Roe.Ave.c).*A{1} +   abs__u.*(Roe.Ave.u.^2+Roe.Ave.v.^2).*A{2}+ abs__u.*Roe.Ave.rho.*Roe.Ave.v.*A{3}+   abs__uc.*(Roe.Ave.H+Roe.Ave.u.*Roe.Ave.c).*A{4})./c22;

   k2 = 2;
   A = Roe.invYdV_G;
   abs__v_c  = EntropyFix( abs(Roe.Ave.v-Roe.Ave.c)).* ProjectedFaceDir{k,k2}  ;
   abs__v    = abs( Roe.Ave.v)                     .* ProjectedFaceDir{k,k2} ; 
   abs__vc   = EntropyFix( abs(Roe.Ave.v+Roe.Ave.c)).* ProjectedFaceDir{k,k2}  ;
   Roe.XLAlpha_G(1,:,:) =  (abs__v_c                               .*A{1} + 2*abs__v                           .*A{3}+                                       abs__vc                               .*A{4})./c22; 
   Roe.XLAlpha_G(3,:,:) =  (abs__v_c.*(Roe.Ave.v-Roe.Ave.c)          .*A{1} + 2*abs__v.*Roe.Ave.v                 .*A{3}+                                       abs__vc.*(Roe.Ave.v+Roe.Ave.c)          .*A{4})./c22; 
   Roe.XLAlpha_G(2,:,:) =  (abs__v_c.*Roe.Ave.u                     .*A{1} + 2*abs__v.*Roe.Ave.u                 .*A{3}+ abs__v.*Roe.Ave.rho          .*A{2}+   abs__vc.*Roe.Ave.u                     .*A{4})./c22;
   Roe.XLAlpha_G(4,:,:) =  (abs__v_c.*(Roe.Ave.H-Roe.Ave.v.*Roe.Ave.c).*A{1} +   abs__v.*(Roe.Ave.u.^2+Roe.Ave.v.^2).*A{3}+ abs__v.*Roe.Ave.rho.*Roe.Ave.u.*A{2}+   abs__vc.*(Roe.Ave.H+Roe.Ave.v.*Roe.Ave.c).*A{4})./c22;

end

function [Roe]=PrepareRoeAverages__given_data_LR(data_L,data_R,Roe)
   global Cp Cv gamma R; 
   sqrho_L = sqrt( data_L.rho );  
   sqrho_R = sqrt( data_R.rho );
   Roe.Ave.rho = sqrho_L.*sqrho_R ;
   Roe.Ave.H   = (data_L.H.*sqrho_L + data_R.H.*sqrho_R) ./ (sqrho_L + sqrho_R) ;
   Roe.Ave.u   = (data_L.u.*sqrho_L + data_R.u.*sqrho_R) ./ (sqrho_L + sqrho_R) ;
   Roe.Ave.v   = (data_L.v.*sqrho_L + data_R.v.*sqrho_R) ./ (sqrho_L + sqrho_R) ;
   Roe.Ave.c   = sqrt( (gamma-1)*(Roe.Ave.H - 0.5*( Roe.Ave.u.^2 + Roe.Ave.v.^2)  ) );
end

function [data_L,data_R] = ReconstructLR_for_RiemannSolverBasedScheme(data_L, data_R, data, method_reconstruction, k)      
    switch  method_reconstruction
        case 'piecewise-constant' %1st order

            if( k ==1 ) 
               data_L.U = data.U(:,1:end-1,:);
               data_R.U = data.U(:,2:end  ,: );
            else
               data_L.U = data.U(:, :,1:end-1);
               data_R.U = data.U(:, :,2:end  );
            end
        otherwise
          errordlg ( ['ReconstructLR: not implemented' method_reconstruction]  ); 
          return
   end
  data_L = U_to_data_LR(data_L);
  data_R = U_to_data_LR(data_R);       
end

function [ AreaFluxA_jhalf, AreaFluxB_jhalf ] = Update_Roe_Flux( data, method_reconstruction  )
    global get_jAjF get_jFjA;
    global get_JFjA_half get_jAJF_half;
    global Sch;
    global Nosqueeze_ProjectedFaceJacob;
    global n_points
    AreaFluxA_jhalf = zeros( 4, n_points(1)-1, n_points(2)); 
    AreaFluxB_jhalf = zeros( 4, n_points(1)  , n_points(2)-1 );
   
    for k = 1:2 %loop two coordiates
       [ Sch{k}.data_L, Sch{k}.data_R] = ReconstructLR_for_RiemannSolverBasedScheme( Sch{k}.data_L, Sch{k}.data_R, data, method_reconstruction, k); 

       % the input vectors for both data_L and data_R is [rho, u, H, P, FluxTerm ] ;
       Sch{k}.Roe =  PrepareRoeAverages__given_data_LR (  Sch{k}.data_L,   Sch{k}.data_R,  Sch{k}.Roe  );
       Sch{k}.Roe =  get_invYdV__given_data_LR( Sch{k}.Roe,  Sch{k}.data_L  ,   Sch{k}.data_R );
       Sch{k}.Roe =  LeftMultiplyBy_RoeAve_X_AbsLambda_PI (Sch{k}.Roe , k );
       
       if k==1 
           AreaFluxA_jhalf(:, :,:) = ( get_JFjA_half(Sch{k}.data_L.F) + get_JFjA_half(Sch{k}.data_R.F)-get_JFjA_half(Sch{k}.Roe.XLAlpha_F) )/2 .* get_jFjA( Nosqueeze_ProjectedFaceJacob{k,1} ) ...
                                     +( get_JFjA_half(Sch{k}.data_L.G) + get_JFjA_half(Sch{k}.data_R.G)-get_JFjA_half(Sch{k}.Roe.XLAlpha_G) )/2 .* get_jFjA( Nosqueeze_ProjectedFaceJacob{k,2} )  ;
       elseif k==2 
           AreaFluxB_jhalf(:, :,: ) = ( get_jAJF_half(Sch{k}.data_L.F) + get_jAJF_half(Sch{k}.data_R.F)- get_jAJF_half(Sch{k}.Roe.XLAlpha_F)  )/2 .* get_jAjF( Nosqueeze_ProjectedFaceJacob{k,1} ) ...
                                     +( get_jAJF_half(Sch{k}.data_L.G) + get_jAJF_half(Sch{k}.data_R.G)- get_jAJF_half(Sch{k}.Roe.XLAlpha_G)  )/2 .* get_jAjF( Nosqueeze_ProjectedFaceJacob{k,2} );
      end

    end % loop two coordiates


end

function [data]= SetInitialConditions(data, BC)
    global Cp Cv gamma R; 
    global n_points; 
    global jAll;
  
    if    contains( BC.type, 'Nozzle'  )

      % Note: here BC.P_R is not used to intialize the last point
      % (j==n_points), however this condition will be prescribed when
      % enforcing the right-side boundary condtions during time-advancement 
      P  = BC.P0_L; 
      T  = BC.T0_L;
      u  = 0;   
      v = 0;
      H  = Cp*T + 0.5 *(u^2 + v^2); 
      rho = P/(T*R) ;
      data.rho( :,:)=rho;  
      data.u(:,:)= u;  
      data.v(:,:)= v; 
      data.H(:,:) = H;
      
      global data_Anlytical_Quasi1D;
      Ny = size(data.rho,2);
      H  = Cp*data_Anlytical_Quasi1D.T + 0.5 *(data_Anlytical_Quasi1D.u.^2); 
      data.rho( :,:)=repmat( data_Anlytical_Quasi1D.rho' , 1 , Ny  );  
      data.u( :,:)=repmat( data_Anlytical_Quasi1D.u' , 1 ,  Ny );  
      data.v(:,:)= 0; 
      data.H(:,:) = repmat( H' , 1 ,  Ny );  
    
    end
    
    
   data = Update_PTEsM_AfterReset_W(data);

end

function [j_Interior, j_BD] = j_Interior_j_BD_for_BoundaryCondtion(k,k2)
    global n_points;
    N = n_points;
    if      (k == 1 &&  k2 ==1)           % bottom-BD
       j_Interior = {2:N(1)-1,  2     };  
       j_BD       = {2:N(1)-1,  1     };     
    elseif  (k == 1 &&  k2==2 )           % top-BD
       j_Interior = {2:N(1)-1,  N(2)-1};   
       j_BD       = {2:N(1)-1,  N(2) };     
    elseif ( k==2  &&  k2==1 )           % left-BD
       j_Interior = {2, 2:N(2)-1};        
       j_BD       = {1, 2:N(2)-1    };     
    elseif ( k==2 &&   k2==2)            %  right-BD
       j_Interior = { N(1)-1, 2:N(2)-1, }; 
       j_BD       =   {N(1), 2:N(2)-1  };     
    end    
end

function [uNormProj,vNormProj,dx_Norm, dx_Parallel] =  Boundary_Projection(k,k2)
   global ProjectedFaceJacob CellVol;
   [j_Interior, j_BD] = j_Interior_j_BD_for_BoundaryCondtion(k,k2) ;
   if( k2 == 1) % bottom or left boundary
      % for normal projection, 3-k is just to flip k 
      dx_Parallel = sqrt( ProjectedFaceJacob{3-k,1}( j_BD{:} ).^2 + ProjectedFaceJacob{3-k,2}( j_BD{:} ).^2)  ;
      uNormProj = ProjectedFaceJacob{3-k,1}( j_BD{:} )./dx_Parallel  ;
      vNormProj = ProjectedFaceJacob{3-k,2}( j_BD{:} )./dx_Parallel  ;
      dx_Norm  = CellVol( j_Interior{:} ) ./ dx_Parallel;  % .* sign(  ProjectedFaceJacob{2,2}( j_BD{:} ) )  ; 
   elseif( k2 == 2) % top or right boundary
      dx_Parallel = sqrt( ProjectedFaceJacob{3-k,1}( j_Interior{:} ).^2 + ProjectedFaceJacob{3-k,2}( j_BD{:} ).^2)  ;
      uNormProj = ProjectedFaceJacob{3-k,1}( j_Interior{:} )./dx_Parallel  ;
      vNormProj = ProjectedFaceJacob{3-k,2}( j_Interior{:} )./dx_Parallel  ;
      dx_Norm  = CellVol( j_Interior{:} ) ./ dx_Parallel;  % .* sign(  ProjectedFaceJacob{2,2}( j_BD{:} ) )  ; 
   end
end

% Specify k(1,2) and k2(1,2) to set a zero-gradient to the "bottom/top/left/right" mesh-boundary
function [data]= BoundaryCondition_ZeroGradient(data,k,k2)
   [j_Interior, j_BD] = j_Interior_j_BD_for_BoundaryCondtion(k,k2);
    data.rho( j_BD{:} ) = data.rho( j_Interior{:} );
    data.u( j_BD{:} ) = data.u( j_Interior{:} );
    data.v( j_BD{:} ) = data.v( j_Interior{:} );
    data.H( j_BD{:} ) = data.H( j_Interior{:} );
end

% Specify k(1,2) and k2(1,2) to set a wall at the "bottom/top/left/right" mesh-boundary
function [data]= BoundaryCondition_Wall(data,dt,k,k2)
    global Cp Cv gamma R; 
    global Nosqueeze_CellVol;
    
    [j_Interior, j_BD]           = j_Interior_j_BD_for_BoundaryCondtion(k,k2);
    [uNormProj,vNormProj,dx_BDNorm , dx_BDPara] =  Boundary_Projection(k,k2);
    
     uNormal_Interior =    data.u(j_Interior{:} ).* uNormProj  + data.v(j_Interior{:})   .* vNormProj ;
     uPara_Interior  = ( data.u(j_Interior{:} ) - uNormal_Interior.* uNormProj );
     vPara_Interior  = ( data.v(j_Interior{:} ) - uNormal_Interior.* vNormProj );
     data.u(j_BD{:}) = uPara_Interior;
     data.v(j_BD{:}) = vPara_Interior;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %set zero wall normal velocity,  zero-gradient to the wall-parallel velocity
     % Force the wall normal-velocity to be zero( a non-zero situation may be caused by inconsistent initialization in the first step) 
     u = data.u( j_BD{:} ) ;
     v = data.v( j_BD{:} ) ;
     uBDNorm= u.* uNormProj  + v.* vNormProj ;
     u = u - uBDNorm .* uNormProj;
     v = v - uBDNorm .* vNormProj ;
     uBDNorm   = 0;
     
     velo_square = u.^2 + v.^2;
     
     rho = data.rho( j_BD{:} );
     H   = data.H  ( j_BD{:} );
     % intermediate values
     T  = ( H-0.5*velo_square )/Cp ; 
     E  = H - R*T;
     P  = rho*R.*T; 
     s  = Cv*log( P./rho.^gamma ); 
     c  = sqrt( gamma * R * T);  

     % The three values representing (negative) time-variation of wave amplitude for each of characertics wave: 
     % L1 = lambda1*(dp/dx-rho*c*du/dx);
     % L2u = lambda2*(c^2*drho/dx-dp/dx)
     % L2v = lambda3*(2*c^2-dv/dx )
     % L3 = lambda4*(dpdx+rho*c*du/dx);
     
     % lambda1-wave goes from the "above"-interior-domain-point to the bottom-wall
     % therefore compute L1 using one-side difference using the interior point
     if( k2 ==1) % bottom or left wall
         L1 = ( uBDNorm-c) .* ( ( data.P(j_Interior{:}) - data.P(j_BD{:})  ) ./dx_BDNorm - rho.* c .* ( uNormal_Interior - uBDNorm )./dx_BDNorm   )  ; 
         % The zero-velocity-wall dictaties the following two conditions
         L3=L1; 
     elseif (k2==2)  % top or right wall
         L3= ( uBDNorm+c ) .* ( ( data.P(j_BD{:}) - data.P(j_Interior{:}) )./dx_BDNorm + rho.*c.* ( uBDNorm -uNormal_Interior) ./dx_BDNorm   )  ; 
         % The zero-velocity-wall dictaties the following two conditions
         L1=L3;         
     end
     L2u = 0;

     
     % The (signed) velocity component parallel to the wall boundary,  here the signed-direction is determined using
     % the dot-Product is between velocity vector [u,v] and the BD direciton vector as [-vNormProj, uNormProj]
     uBDPara = sqrt(velo_square) .*  sign( u.*(vNormProj) + v.* (-uNormProj)   ) ; 
     % Now update equation of continiuty and energy ( momentum equation is not needed)
     d1 = (L2u+ (L3+L1)/2 )./c.^2; 
     d2 = (L3+L1)/2; 
     d3 = (L3-L1)./(2*rho.*c);

     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     rho_new  = rho    - dt* d1                 ;
     rhoE_new = rho.*E - dt* (0.5*velo_square.*d1+  d2/(gamma-1) + rho.*uBDNorm.*d3 )    ;
     
     rhouBDPara_new = rho.*uBDPara -dt*( uBDPara.*d1 )                                 ;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     


     %%%%%% caculation using interior scheme  %%%%%%
     %data.U(:, j_BD{:} ) = data.U(:, j_BD{:} )- dt* ( data.AreaFluxA_jhalf(:,j_BD{:})-data.AreaFluxA_jhalf(:,j_BD{1}-1, j_BD{2})    +data.AreaFluxB_jhalf(:,j_BD{:})- data.AreaFluxB_jhalf(:,j_BD{1},j_BD{1}-1)      ) ./ Nosqueeze_CellVol(:, j_BD{:} );     
     if( k == 1)  % bottom/top wall
        out =  - dt* (    data.AreaFluxA_jhalf(:,j_BD{:})-data.AreaFluxA_jhalf(:,j_BD{1}-1, j_BD{2})          ) ./ Nosqueeze_CellVol(:, j_BD{:} );      
     elseif (k==2) % left/right wall 
        out =  - dt* (    data.AreaFluxB_jhalf(:,j_BD{:}) - data.AreaFluxB_jhalf(:,j_BD{1},j_BD{2}-1) ) ./ Nosqueeze_CellVol(:, j_BD{:} );      
     end     

     z = size(out);
     rho_new  = rho_new +   reshape( out(1,:) , [ z(2:end) 1] ) ; 
     rhoE_new = rhoE_new +  reshape( out(4,:),  [ z(2:end) 1] ) ; 
     rhouBDPara_new  = rhouBDPara_new +  reshape( out(2,:),  [ z(2:end) 1] ).* vNormProj +  reshape( out(3,:),  [ z(2:end) 1] ).* (-uNormProj);

     % the newly updated condtions 
     uBDPara_new = rhouBDPara_new./rho_new;
     
     %uBDPara_new = uBDPara;
     
     E_new   =  rhoE_new./rho_new; 
     T_new   = (E_new  -0.5*uBDPara_new.^2 ) /Cv; %
     H_new   = E_new + R*T_new;
     
     data.rho(j_BD{:}) = rho_new;
     data.H(j_BD{:})   = H_new; 

     data.u(j_BD{:})= uBDPara_new.*vNormProj;   % 2D projection of BD paralle velocity
     data.v(j_BD{:})= uBDPara_new.*(-uNormProj);% 2D projection of BD paralle velocity

end

function [data]= BoundaryCondition_Nozzle_LeftSubsonicInlet_RightOutlet(data,BC)
    global Cp Cv gamma R; 
    global n_points; 
    global ProjectedFaceJacob ;
    
    N = n_points;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % !!!!!! left and right booundary !!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  For 1D euler equation there are three (sorted) eigenvalues:  
    %                  u-c, u , u+c
    %     which corresponds to three characteristical variables: 
    %   1:  characteristic_var_first = c/(gamma-1)-u/2 
    %   2:  characteristic_var_second= s                  (i.e. entropy ) 
    %   3:  characteristic_var_third = c/(gamma-1)+u/2
    %
    % The following set boundary condition for the leftmost cell 1 which is always assumed to
    %    be subsonic, therefore there are  two waves going into domain, one coming out the 
    %    domain from the neighbouring cell 2,  we choose to specify two variables of P0_L and T0_L: 
    u = 2 * data.u(2, : ) - data.u(3,: ) ;
    v = 2 * data.v(2, : ) - data.v(3,: ) ;
    velo = sqrt(u.^2 + v.^2) ;
    H = 2 * data.H(2, : ) - data.H(3,: ) ;
    T = (H-0.5*velo.^2 )/Cp ; 
    c = sqrt( gamma * R * T ) ;
    characteristic_var_first = c /(gamma-1) - velo/2; 
    
    % solve out [rho,u,H] from three given values of (1) P0_L (2) T0_L and (3) characteristic_var_first 
    % Note: the relation for "total" temperature and pressure 
    %        T0/T = 1+ 0.5*(gamma-1)*M^2  = 1 + 0.5* (gamma-1)*(u/c)^2
    %        P0/P = (T0/T )^(1/gamma -1 ) 
    %          Here is a 2nd order polynomial equation for the speed of sound c (let x1 denotes  characteristic_var_first)
    %           [1+2/(gamma-1)]* c^2 - [4*x1] *c +  [ 2(gamma-1) x1^2 - T0*gamma*R  ] = 0  
    a1 =1+2/(gamma-1);  a2 = -4 * characteristic_var_first; a3 = 2*(gamma-1)*characteristic_var_first.^2 - BC.T0_L*gamma*R; 
    c  = (-a2 + sqrt ( a2.^2 - 4*a1.*a3) )./(2*a1);  % speed of sound
    T  = (c.^2)/(gamma*R);
    P  = BC.P0_L ./ ( (BC.T0_L./T).^(gamma/(gamma -1 ) ) );
    rho= P./(R*T );
    velo  = 2* ( sqrt(gamma*R*T) /(gamma-1)  -  characteristic_var_first ) ; 
    data.rho(1,:) = rho ; 
    data.u(1,:)   =  velo .* ProjectedFaceJacob{1,1}(1,:)  ./ sqrt( ProjectedFaceJacob{1,1}(1,:).^2 + ProjectedFaceJacob{1,2}(1,:).^2) ; 
    data.v(1,:)   =  velo .* ProjectedFaceJacob{1,2}(1,:)  ./ sqrt( ProjectedFaceJacob{1,1}(1,:).^2 + ProjectedFaceJacob{1,2}(1,:).^2) ; 
    data.H(1,:)   = Cp*T + 0.5*velo.^2; 
    %
    %  The following code set the boundary condition on the rightmost cell N
    %    which can be either subsonic or supersonic, 
    %  supersonic: nothing is needed at N
    %  subsonic  : one condition is required to be prescibed at N, we choose it to be P_R
    rho =2* data.rho(end-1,:) - data.rho(end-2,:);
    u   =2* data.u(end-1,:)   - data.u(end-2,:) ;
    v   =2* data.v(end-1,:)   - data.v(end-2,:) ;
    H   =2* data.H(end-1,:)   - data.H(end-2,:);
    velo = sqrt(u.^2 + v.^2) ;

    T = (H-0.5*velo.^2 )/Cp ; 
    P = rho*R.*T; 
    s = Cv*log( P./rho.^gamma ); 
    c = sqrt( gamma * R *T ) ;

    characteristic_Var_second =  s ;    
    characteristic_Var_third  =  c/(gamma-1)+velo/2  ;
    % solve [rho,u,H] from three given values of (1) P_R (2) characteristic_Var_second and (3) characteristic_Var_third 

    P   = BC.P_R;
    rho = power( P ./ ( exp( characteristic_Var_second/Cv ) )  , 1./gamma  ) ; 
    T   = P ./( rho * R ) ; 
    c   = sqrt(gamma*R*T);
    velo   = 2* ( characteristic_Var_third - c/( gamma-1) ); 

    %[uNormProj,vNormProj,dx_BDNorm ] =  Boundary_Projection(1,2);
    
    data.rho(end,:)  =  rho;
    data.u(end,:)    =  velo .* ProjectedFaceJacob{1,1}(end-1,:)  ./ sqrt( ProjectedFaceJacob{1,1}(end-1,:).^2 + ProjectedFaceJacob{1,2}(end-1,:).^2) ; 
    data.v(end,:)    =  velo .* ProjectedFaceJacob{1,2}(end-1,:)  ./ sqrt( ProjectedFaceJacob{1,1}(end-1,:).^2 + ProjectedFaceJacob{1,2}(end-1,:).^2) ; 
    data.H(end,:)    = Cp*T+0.5*velo.^2;
end

function [data]= BoundaryCondition_SupersonicInlet(data,BC,k,k2)
    global Cp Cv gamma R;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left, Supersonic inlet
    if( BC.M_L < 1) 
       errordlg( 'BC.M_L < 1 for supersonic inlet');
       return
    end       
    
    velo = sqrt( gamma*R*BC.T_L)*BC.M_L; 
    H  = Cp*BC.T_L + 0.5 *(velo^2);  
    rho = BC.P_L/( R*BC.T_L ) ; 
 
    [j_Interior, j_BD] = j_Interior_j_BD_for_BoundaryCondtion(k,k2) ;
    [uNormProj,vNormProj,dx_BDNorm, dx_Parallel ] =  Boundary_Projection(k,k2);
    data.rho( j_Interior{:} )=rho;   
    data.u(j_Interior{:})= velo*uNormProj;
    data.v(j_Interior{:})= velo*vNormProj;
    data.H(j_Interior{:}) = H;
end

function Prepare_and_AllocateMemoryForComputionalSpeedUp()
   global get_jAjF get_jFjA get_jF1jA get_jAjF1 get_jF1jo     get_jFjo       get_jojF       get_jojF1    get_jojo;
   global get_JFjA_half get_jAJF_half    get_Jojo_half get_joJo_half  get_J_1jo_half get_joJ_1_half;
   global n_points;
   global data;
   global data_nhalf;
   global Sch;
   get_jAjF  = @(UFG) UFG(:,:,1:end-1 );
   get_jAjF1 = @(UFG) UFG(:,:,2:end );
   
   get_jFjA  = @(UFG) UFG(:,1:end-1,:);
   get_jF1jA = @(UFG) UFG(:,2:end,  :);
   
   
   get_jF1jo = @(UFG) UFG(:,2:end,  2:end-1 );
   get_jFjo  = @(UFG) UFG(:,1:end-1,2:end-1 );
   get_jojF  = @(UFG) UFG(:,2:end-1,1:end-1 );
   get_jojF1 = @(UFG) UFG(:,2:end-1,2:end );
   
   get_jojo  = @(UFG) UFG(:,2:end-1,2:end-1 );
  
   
   get_JFjA_half  = @(UFG) UFG(:,:,:);
   get_jAJF_half  = @(UFG) UFG(:,:,: );
   
   get_Jojo_half  = @(UFG) UFG(:,2:end,2:end-1);
   get_joJo_half  = @(UFG) UFG(:,2:end-1,2:end);
   
   get_J_1jo_half = @(UFG) UFG(:,1:end-1,2:end-1 );
   get_joJ_1_half = @(UFG) UFG(:,2:end-1,1:end-1 );  


   
   
   data.u = zeros( n_points(1), n_points(2) );
   data.v = zeros( n_points(1), n_points(2) );
   data.P = zeros( n_points(1), n_points(2) );
   data.T = zeros( n_points(1), n_points(2) );
   data.H = zeros( n_points(1), n_points(2) );
   data.rho = zeros( n_points(1), n_points(2) );
   data.E  = zeros( n_points(1), n_points(2) );
   data.s  = zeros( n_points(1), n_points(2) );
   data.U  = zeros( 4, n_points(1), n_points(2) );
   data.F  = zeros( 4, n_points(1), n_points(2) );
   data.G  = zeros( 4, n_points(1), n_points(2) );
   data.F_jhalf  = zeros(4, n_points(1)-1, n_points(2)   );
   data.G_jhalf  = zeros(4, n_points(1)  , n_points(2)-1 );
   data.AreaFluxA_jhalf  = zeros( 4, n_points(1)-1, n_points(2)   );
   data.AreaFluxB_jhalf  = zeros( 4, n_points(1)  , n_points(2)-1 );

   
   data_nhalf = data;

   
   for k = 1: 2
      if  k==1
         N = { n_points(1)-1, n_points(2)  };
      elseif k==2
         N = { n_points(1), n_points(2)-1  };
      end
      
      Sch{k}.data_L.u = zeros( N{:} );
      Sch{k}.data_L.v = zeros( N{:} );
      Sch{k}.data_L.P = zeros( N{:} );
      Sch{k}.data_L.T = zeros( N{:} );
      Sch{k}.data_L.H = zeros( N{:} );
      Sch{k}.data_L.rho = zeros( N{:} );
      Sch{k}.data_L.E  = zeros( N{:} );
      %Sch{k}.data_L.s  = zeros( N{:} );
      Sch{k}.data_L.U  = zeros( 4, N{:} );
      Sch{k}.data_L.F  = zeros( 4, N{:} );
      Sch{k}.data_L.G  = zeros( 4, N{:} );


      Sch{k}.data_R.u = zeros( N{:} );
      Sch{k}.data_R.v = zeros( N{:} );
      Sch{k}.data_R.P = zeros( N{:} );
      Sch{k}.data_R.T = zeros( N{:});
      Sch{k}.data_R.H = zeros( N{:} );
      Sch{k}.data_R.rho = zeros( N{:} );
      Sch{k}.data_R.E  = zeros( N{:} );
      %Sch{k}.data_R.s  = zeros( N{:} );
      Sch{k}.data_R.U  = zeros( 4, N{:} );
      Sch{k}.data_R.F  = zeros( 4, N{:});
      Sch{k}.data_R.G  = zeros( 4, N{:});

      Sch{k}.Roe.XLAlpha_F   = zeros(4, N{:} );
      Sch{k}.Roe.XLAlpha_G   = zeros(4, N{:}  );
      for m =1 : 4
         Sch{k}.Roe.invYdV_F{m}    = zeros( N{:} );
         Sch{k}.Roe.invYdV_G{m}    = zeros( N{:} );
      end
      Sch{k}.Roe.Ave.rho     =  zeros( N{:} );
      Sch{k}.Roe.Ave.H       =  zeros( N{:} );
      Sch{k}.Roe.Ave.u       =  zeros( N{:} );
      Sch{k}.Roe.Ave.v       =  zeros( N{:} );
      Sch{k}.Roe.Ave.c       =  zeros( N{:} );
   end

end

