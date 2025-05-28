% -------------------------------------------------
%  2D nozzle data generation using pressure splits 
% -------------------------------------------------
clear; clc; close all;


%% 1) Desired number of valid samples & parameter ranges
nSamplesDesired = 10;  % total # of samples

nSubsonic   = floor(nSamplesDesired / 3);  
nShock      = floor(nSamplesDesired / 3);  
nSupersonic = nSamplesDesired - nSubsonic - nShock;  

% Stagnation pressure
P0_min = 1.5e6;  % [Pa]
P0_max = 2e6;

% T0 range
T0_min = 280;   % [K]
T0_max = 400;   % [K]

% Three back-pressure bands
% Subsonic
PR_min_sub   = 1.85e6;
PR_max_sub   = 1.99e6;

% Shock
PR_min_shock = 8e5;
PR_max_shock = 1.85e6;

% Supersonic
PR_min_super = 1e5;
PR_max_super = 8e5;

% Geometry ranges
% rIn_min = 0.65;   rIn_max = 1.0;
% rTh_min = 0.3;    rTh_max = 0.5;
% rOut_min= 0.65;   rOut_max= 0.85;
% xTh_min = -1.5;   xTh_max = 1.5;
rIn_min = 0.7;   rIn_max = 1;
rTh_min = 0.3;    rTh_max = 0.5;
rOut_min= 0.7;   rOut_max= 1;
xTh_min = -1.5;   xTh_max = 1.5;

%% 2) LHS for Subsonic band
nLarge_sub   = 15 * nSubsonic;
lhsPointsSub = lhsdesign(nLarge_sub, 7);  % columns: {T0,P0, PR, rIn, rTh, rOut, xTh}

T0_sub_raw   = T0_min + (T0_max - T0_min) .* lhsPointsSub(:,1);
P0_sub_raw   = P0_min + (P0_max - P0_min).*lhsPointsSub(:,2);
PR_sub_raw   = PR_min_sub + (PR_max_sub - PR_min_sub) .* lhsPointsSub(:,3);
rIn_sub_raw  = rIn_min  + (rIn_max  - rIn_min)  .* lhsPointsSub(:,4);
rTh_sub_raw  = rTh_min  + (rTh_max  - rTh_min)  .* lhsPointsSub(:,5);
rOut_sub_raw = rOut_min + (rOut_max - rOut_min).* lhsPointsSub(:,6);
xTh_sub_raw  = xTh_min  + (xTh_max  - xTh_min)  .* lhsPointsSub(:,7);

validMask_sub = (PR_sub_raw < P0_sub_raw & P0_sub_raw - PR_sub_raw > 50000);
T0_sub_raw   = T0_sub_raw(validMask_sub);
P0_sub_raw = P0_sub_raw(validMask_sub);
PR_sub_raw   = PR_sub_raw(validMask_sub);
rIn_sub_raw  = rIn_sub_raw(validMask_sub);
rTh_sub_raw  = rTh_sub_raw(validMask_sub);
rOut_sub_raw = rOut_sub_raw(validMask_sub);
xTh_sub_raw  = xTh_sub_raw(validMask_sub);

if numel(T0_sub_raw) < nSubsonic
    error('Not enough valid subsonic combos. Adjust subsonic band or sampling.');
end

T0_sub   = T0_sub_raw(1:nSubsonic);
P0_sub = P0_sub_raw(1:nSubsonic);
PR_sub   = PR_sub_raw(1:nSubsonic);
rIn_sub  = rIn_sub_raw(1:nSubsonic);
rTh_sub  = rTh_sub_raw(1:nSubsonic);
rOut_sub = rOut_sub_raw(1:nSubsonic);
xTh_sub  = xTh_sub_raw(1:nSubsonic);

%% 3) LHS for Shock band
nLarge_shock   = 15 * nShock;
lhsPointsShock = lhsdesign(nLarge_shock, 7);

T0_sh_raw   = T0_min + (T0_max - T0_min) .* lhsPointsShock(:,1);
PR_sh_raw   = PR_min_shock + (PR_max_shock - PR_min_shock).* lhsPointsShock(:,2);
P0_sh_raw   = P0_min + (P0_max - P0_min).* lhsPointsShock(:,3);
rIn_sh_raw  = rIn_min  + (rIn_max  - rIn_min)  .* lhsPointsShock(:,4);
rTh_sh_raw  = rTh_min  + (rTh_max  - rTh_min)  .* lhsPointsShock(:,5);
rOut_sh_raw = rOut_min + (rOut_max - rOut_min).* lhsPointsShock(:,6);
xTh_sh_raw  = xTh_min  + (xTh_max  - xTh_min)  .* lhsPointsShock(:,7);

validMask_sh = (PR_sh_raw < P0_sh_raw & P0_sh_raw - PR_sh_raw > 50000);
T0_sh_raw   = T0_sh_raw(validMask_sh);
PR_sh_raw   = PR_sh_raw(validMask_sh);
P0_sh_raw   = P0_sh_raw(validMask_sh);
rIn_sh_raw  = rIn_sh_raw(validMask_sh);
rTh_sh_raw  = rTh_sh_raw(validMask_sh);
rOut_sh_raw = rOut_sh_raw(validMask_sh);
xTh_sh_raw  = xTh_sh_raw(validMask_sh);

if numel(T0_sh_raw) < nShock
    error('Not enough valid shock combos. Adjust shock band or sampling.');
end

T0_sh   = T0_sh_raw(1:nShock);
PR_sh   = PR_sh_raw(1:nShock);
P0_sh   = P0_sh_raw(1:nShock);
rIn_sh  = rIn_sh_raw(1:nShock);
rTh_sh  = rTh_sh_raw(1:nShock);
rOut_sh = rOut_sh_raw(1:nShock);
xTh_sh  = xTh_sh_raw(1:nShock);

%% 4) LHS for Supersonic band
nLarge_super   = 15 * nSupersonic;
lhsPointsSuper = lhsdesign(nLarge_super, 7);

T0_sup_raw   = T0_min + (T0_max - T0_min) .* lhsPointsSuper(:,1);
PR_sup_raw   = PR_min_super + (PR_max_super - PR_min_super).* lhsPointsSuper(:,2);
P0_sup_raw   = P0_min + (P0_max - P0_min).* lhsPointsSuper(:,3);
rIn_sup_raw  = rIn_min  + (rIn_max  - rIn_min)  .* lhsPointsSuper(:,4);
rTh_sup_raw  = rTh_min  + (rTh_max  - rTh_min)  .* lhsPointsSuper(:,5);
rOut_sup_raw = rOut_min + (rOut_max - rOut_min).* lhsPointsSuper(:,6);
xTh_sup_raw  = xTh_min  + (xTh_max  - xTh_min)  .* lhsPointsSuper(:,7);

validMask_sup = (PR_sup_raw < P0_sup_raw & P0_sup_raw - PR_sup_raw > 50000);
T0_sup_raw   = T0_sup_raw(validMask_sup);
PR_sup_raw   = PR_sup_raw(validMask_sup);
P0_sup_raw   = P0_sup_raw(validMask_sup);
rIn_sup_raw  = rIn_sup_raw(validMask_sup);
rTh_sup_raw  = rTh_sup_raw(validMask_sup);
rOut_sup_raw = rOut_sup_raw(validMask_sup);
xTh_sup_raw  = xTh_sup_raw(validMask_sup);

if numel(T0_sup_raw) < nSupersonic
    error('Not enough valid supersonic combos. Adjust super band or sampling.');
end

T0_sup   = T0_sup_raw(1:nSupersonic);
PR_sup   = PR_sup_raw(1:nSupersonic);
P0_sup   = P0_sup_raw(1:nSupersonic);
rIn_sup  = rIn_sup_raw(1:nSupersonic);
rTh_sup  = rTh_sup_raw(1:nSupersonic);
rOut_sup = rOut_sup_raw(1:nSupersonic);
xTh_sup  = xTh_sup_raw(1:nSupersonic);

%% 5) Merge all
T0_all   = [T0_sub;   T0_sh;   T0_sup];
PR_all   = [PR_sub;   PR_sh;   PR_sup];
P0_all   = [P0_sub;   P0_sh;   P0_sup];
rIn_all  = [rIn_sub;  rIn_sh;  rIn_sup];
rTh_all  = [rTh_sub;  rTh_sh;  rTh_sup];
rOut_all = [rOut_sub; rOut_sh; rOut_sup];
xTh_all  = [xTh_sub;  xTh_sh;  xTh_sup];

nTotal = length(T0_all);
if nTotal ~= nSamplesDesired
    error('We ended up with %d total samples, expected %d.', nTotal, nSamplesDesired);
end

% Shuffle so they aren't in band-order
idxShuffle = randperm(nTotal);
T0_all   = T0_all(idxShuffle);
PR_all   = PR_all(idxShuffle);
P0_all   = P0_all(idxShuffle);
rIn_all  = rIn_all(idxShuffle);
rTh_all  = rTh_all(idxShuffle);
rOut_all = rOut_all(idxShuffle);
xTh_all  = xTh_all(idxShuffle);


%% 6) Domain
DomainSize_in_x = 5;  
n_points_x        = 200;  
n_points_y = 50;
x_1D = linspace(-DomainSize_in_x/2, DomainSize_in_x/2, n_points_x)';
%% 7) Counters
shockInsideCount = 0;
shockAtExitCount = 0;
subsonicCount    = 0;


%% 8) Preallocate BC data + solution arrays
% Store BC columns: [T0_L, P0_L, P_R, x_throat, r_inlet, r_throat, r_out]
bcData = zeros(nSamplesDesired, 7);

% Each solution array is [nSamplesDesired × n_points]
X_all            = zeros(nSamplesDesired, n_points_x, n_points_y );
Y_all            = zeros(nSamplesDesired, n_points_x, n_points_y);
rho_all          = zeros(nSamplesDesired, n_points_x, n_points_y);
T_all_solutions  = zeros(nSamplesDesired, n_points_x, n_points_y);
P_all_solutions  = zeros(nSamplesDesired, n_points_x, n_points_y);
u_all            = zeros(nSamplesDesired, n_points_x, n_points_y);
v_all            = zeros(nSamplesDesired, n_points_x, n_points_y);
%% 9) Main loop over samples
for i = 1 : nSamplesDesired
    % Build BC struct
    T0_current = T0_all(i);
    P0_current = P0_all(i);
    PR_current = PR_all(i);
    bcData(i,:) = [ T0_current, P0_current, PR_current, ...
                    xTh_all(i), rIn_all(i), rTh_all(i), rOut_all(i) ];
    fprintf('\nSample %d:\n', i);
    fprintf('  T0_L     = %.2f K\n',  T0_current);
    fprintf('  P0_L     = %.2f Pa\n', P0_current);
    fprintf('  P_R      = %.2f Pa\n', PR_current);
    fprintf('  x_throat = %.4f\n',    xTh_all(i));
    fprintf('  r_inlet  = %.4f\n',    rIn_all(i));
    fprintf('  r_throat = %.4f\n',    rTh_all(i));
    fprintf('  r_out    = %.4f\n',    rOut_all(i));


    % Build nozzle geometry
    r_nozzle = zeros(size(x_1D));
    [convergingFun, divergingFun] = nozzleGeom(...
        -DomainSize_in_x/2, xTh_all(i), DomainSize_in_x/2, ...
         rIn_all(i),         rTh_all(i), rOut_all(i)        );

    for j = 1:n_points_x
        xx = x_1D(j);
        if xx <= xTh_all(i)
            r_nozzle(j) = convergingFun(xx);
        else
            r_nozzle(j) = divergingFun(xx);
        end
    end

    % Quick checks
    idx_converging = (x_1D <= xTh_all(i));
    idx_diverging  = (x_1D > xTh_all(i));
    if any(diff(r_nozzle(idx_converging)) > 0)
        error('Converging section not strictly convergent (sample %d).', i);
    end
    if any(diff(r_nozzle(idx_diverging)) < 0)
        error('Diverging section not strictly divergent (sample %d).', i);
    end
    

    % Solve the flow
    [data, x, y, solverBranch] = Euler2D_FunctionV3(P0_current, PR_current, T0_current, r_nozzle);

    % Count solution type
    switch solverBranch
        case 'supersonicShock'
            shockInsideCount = shockInsideCount + 1;
        case 'supersonic'
            shockAtExitCount = shockAtExitCount + 1;
        case 'subsonic'
            subsonicCount    = subsonicCount + 1;
        otherwise
            warning('Unknown solver branch: %s', solverBranch);
    end
    X_all(i,:, :)           = x;
    Y_all(i,:,:)            = y;
    rho_all(i,:,:)          = data.rho;
    T_all_solutions(i,:,:)  = data.T;
    P_all_solutions(i,:,:)  = data.P;
    u_all(i,:,:)            = data.u;
    v_all(i,:,:)            = data.v;
end

   

%% 10) Print distribution of 1D quasi solutions
fprintf('\n\n---------------------------------------\n');
fprintf('Number of shock (internal) solutions : %d\n', shockInsideCount);
fprintf('Number of fully supersonic solutions : %d\n', shockAtExitCount);
fprintf('Number of fully subsonic solutions   : %d\n', subsonicCount);
fprintf('Data generation complete!\n');

%% 12) Write everything to a single HDF5 file
h5File1 = 'C:\Thesis\Python FNO\2DTest_x_y.h5';
h5File2 = 'C:\Thesis\Python FNO\2DTestT_P.h5';
h5File3 = 'C:\Thesis\Python FNO\2DTestu_v.h5';
if exist(h5File1, 'file')
    delete(h5File1);  % Overwrite if it already exists
end


if exist(h5File2, 'file')
    delete(h5File2);  % Overwrite if it already exists
end


if exist(h5File3, 'file')
    delete(h5File3);  % Overwrite if it already exists
end

% bcData has shape [nSamplesDesired, 7].
% Store it under the dataset name "/BC".
h5create(h5File1, '/BC', size(bcData'));
h5write(h5File1, '/BC', bcData');

% Each solution array is [nSamplesDesired × n_points].
h5create(h5File1, '/X', size(X_all));
h5write(h5File1, '/X', X_all);

h5create(h5File1, '/Y', size(Y_all));
h5write(h5File1, '/Y', Y_all);

h5create(h5File2, '/T', size(T_all_solutions));
h5write(h5File2, '/T', T_all_solutions);

h5create(h5File2, '/P', size(P_all_solutions));
h5write(h5File2, '/P', P_all_solutions);

h5create(h5File3, '/u', size(u_all));
h5write(h5File3, '/u', u_all);

h5create(h5File3, '/v', size(v_all));
h5write(h5File3, '/v', v_all);

