%---------------------------------------------------------
% 1D quasi nozzle data generation using pressure splits,
% storing all BC and Output data in a single HDF5 file.
%---------------------------------------------------------

clear; clc; close all;

%% 1) Global constants
global Cp Cv gamma R;
Cp    = 1000;
gamma = 1.4;
Cv    = Cp / gamma;
R     = Cp - Cv;

%% 2) Desired number of valid samples & parameter ranges
nSamplesDesired = 200;  % total # of samples

% ~1/3 in subsonic band, ~1/3 shock band, ~1/3 supersonic band:
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
rIn_min = 0.7;   rIn_max = 1;
rTh_min = 0.3;    rTh_max = 0.5;
rOut_min= 0.7;   rOut_max= 1;
xTh_min = -1.5;   xTh_max = 1.5;

%% 3) LHS for Subsonic band
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

%% 4) LHS for Shock band
nLarge_shock   = 10 * nShock;
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

%% 5) LHS for Supersonic band
nLarge_super   = 10 * nSupersonic;
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

%% 6) Merge all
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

% (Optional) Shuffle so they aren't in band-order
idxShuffle = randperm(nTotal);
T0_all   = T0_all(idxShuffle);
PR_all   = PR_all(idxShuffle);
P0_all   = P0_all(idxShuffle);
rIn_all  = rIn_all(idxShuffle);
rTh_all  = rTh_all(idxShuffle);
rOut_all = rOut_all(idxShuffle);
xTh_all  = xTh_all(idxShuffle);

%% 7) Domain setup
DomainSize_in_x = 5;
n_points        = 128;
x_1D = linspace(-DomainSize_in_x/2, DomainSize_in_x/2, n_points)';

%% 8) Initialize counters
shockInsideCount = 0;
shockAtExitCount = 0;
subsonicCount    = 0;

%% 9) Preallocate BC data + solution arrays
% Store BC columns: [T0_L, P0_L, P_R, x_throat, r_inlet, r_throat, r_out]
bcData = zeros(nSamplesDesired, 7);

% Each solution array is [nSamplesDesired × n_points]
X_all            = zeros(nSamplesDesired, n_points);
Y_all            = zeros(nSamplesDesired, n_points);
rho_all          = zeros(nSamplesDesired, n_points);
T_all_solutions  = zeros(nSamplesDesired, n_points);
P_all_solutions  = zeros(nSamplesDesired, n_points);
u_all            = zeros(nSamplesDesired, n_points);

%% 10) Main loop over samples
for i = 1 : nSamplesDesired
    % Boundary conditions
    T0_current = T0_all(i);
    P0_current = P0_all(i);
    PR_current = PR_all(i);

    % Populate bcData row (i):
    % [T0_L, P0_L, P_R, x_throat, r_inlet, r_throat, r_out]
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
    [convergingFun, divergingFun] = nozzleGeom( ...
        -DomainSize_in_x/2, xTh_all(i), DomainSize_in_x/2, ...
         rIn_all(i),        rTh_all(i), rOut_all(i));

    for j = 1:n_points
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
    area = pi * (r_nozzle.^2);
    BC.type = 'VenturiNozzleFlow';
    BC.T0_L = T0_current;
    BC.P0_L = P0_current;
    BC.P_R  = PR_current;

    [data_Anlytical_Quasi1D, solverBranch] = ExactSteadyStateSolution_Quasi1DNozzleFlow_Max(n_points, r_nozzle, BC);

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

    % Clean up any NaNs in the solution
    fieldsSol = fieldnames(data_Anlytical_Quasi1D);
    for f = 1 : numel(fieldsSol)
        vec = data_Anlytical_Quasi1D.(fieldsSol{f});
        if isvector(vec) && length(vec) == n_points
            nanIdx = isnan(vec);
            if any(nanIdx)
                xVals      = (1:length(vec))';
                vec(nanIdx) = interp1(xVals(~nanIdx), vec(~nanIdx), xVals(nanIdx), 'linear','extrap');
                data_Anlytical_Quasi1D.(fieldsSol{f}) = vec;
            end
        end
    end

    % Store outputs in our big arrays
    X_all(i,:)            = x_1D(:)';
    Y_all(i,:)            = r_nozzle(:)';
    rho_all(i,:)          = data_Anlytical_Quasi1D.rho(:)';
    T_all_solutions(i,:)  = data_Anlytical_Quasi1D.T(:)';
    P_all_solutions(i,:)  = data_Anlytical_Quasi1D.P(:)';
    u_all(i,:)            = data_Anlytical_Quasi1D.u(:)';
end

%% 11) Print distribution of solutions
fprintf('\n\n---------------------------------------\n');
fprintf('Number of shock (internal) solutions : %d\n', shockInsideCount);
fprintf('Number of fully supersonic solutions : %d\n', shockAtExitCount);
fprintf('Number of fully subsonic solutions   : %d\n', subsonicCount);
fprintf('Data generation complete!\n');

%% 12) Write everything to a single HDF5 file
h5File = 'C:\Thesis\Python FNO\1D15000V5.h5';
if exist(h5File, 'file')
    delete(h5File);  % Overwrite if it already exists
end

% bcData has shape [nSamplesDesired, 7].
% Store it under the dataset name "/BC".
h5create(h5File, '/BC', size(bcData'));
h5write(h5File, '/BC', bcData');

% Each solution array is [nSamplesDesired × n_points].
h5create(h5File, '/X', size(X_all'));
h5write(h5File, '/X', X_all');

h5create(h5File, '/Y', size(Y_all'));
h5write(h5File, '/Y', Y_all');

h5create(h5File, '/rho', size(rho_all'));
h5write(h5File, '/rho', rho_all');

h5create(h5File, '/T', size(T_all_solutions'));
h5write(h5File, '/T', T_all_solutions');

h5create(h5File, '/P', size(P_all_solutions'));
h5write(h5File, '/P', P_all_solutions');

h5create(h5File, '/u', size(u_all'));
h5write(h5File, '/u', u_all');

fprintf('All data written to: %s\n', h5File);