% ======================================
%  Piecewise Converging-Diverging Nozzle 
% ======================================


function [convergingFun, divergingFun] = nozzleGeom(x_inlet, x_throat, x_outlet, R_inlet, R_throat, R_outlet)
    
    
    % Define two cubic polynomials:
    %   r1(x) on [x_inlet, x_throat], 
    %   r2(x) on [x_throat, x_outlet].
    
    %-------------------------
    % 2) Converging Section Coeffs: r1(x) = a0 + a1*x + a2*x^2 + a3*x^3
    %    Domain: x_inlet <= x <= x_throat   i.e.  -2.5 <= x <= 0
    %    Conditions:
    %      r1(x_inlet)  = R_inlet
    %      r1'(x_inlet) = 0
    %      r1(x_throat) = R_throat
    %      r1'(x_throat)= 0
    %-------------------------
    
    syms a0 a1 a2 a3 real
    
    % r1(x) and its derivative:
    r1  = @(x) a0 + a1*x + a2*x^2 + a3*x^3;
    dr1 = @(x) a1 + 2*a2*x + 3*a3*x^2;
    
    eq1 = r1(x_inlet)  == R_inlet;    % r1(-2.5) = R_inlet
    eq2 = dr1(x_inlet) == 0;          % r1'(-2.5) = 0
    eq3 = r1(x_throat) == R_throat;   % r1(0) = R_throat
    eq4 = dr1(x_throat)== 0;          % r1'(0) = 0
    
    sol1 = solve([eq1, eq2, eq3, eq4],[a0, a1, a2, a3], 'Real', true);
    
    a0_1 = double(sol1.a0);
    a1_1 = double(sol1.a1);
    a2_1 = double(sol1.a2);
    a3_1 = double(sol1.a3);
    
    % For clarity, define a function handle for the converging section:
    convergingFun = @(x) a0_1 + a1_1*x + a2_1*x.^2 + a3_1*x.^3;
    
    %-------------------------
    % 3) Diverging Section Coeffs: r2(x) = b0 + b1*x + b2*x^2 + b3*x^3
    %    Domain: x_throat <= x <= x_outlet   i.e.  0 <= x <= 2.5
    %    Conditions:
    %      r2(x_throat) = R_throat
    %      r2'(x_throat)= 0
    %      r2(x_outlet) = R_outlet
    %      r2'(x_outlet)= 0
    %-------------------------
    
    syms b0 b1 b2 b3 real
    
    r2  = @(x) b0 + b1*x + b2*x^2 + b3*x^3;
    dr2 = @(x) b1 + 2*b2*x + 3*b3*x^2;
    
    eq5  = r2(x_throat)  == R_throat;   % r2(0) = R_throat
    eq6  = dr2(x_throat) == 0;          % r2'(0) = 0
    eq7  = r2(x_outlet)  == R_outlet;   % r2(2.5)= R_outlet
    eq8  = dr2(x_outlet) == 0;          % r2'(2.5)=0
    
    sol2 = solve([eq5, eq6, eq7, eq8],[b0,b1,b2,b3], 'Real', true);
    
    b0_2 = double(sol2.b0);
    b1_2 = double(sol2.b1);
    b2_2 = double(sol2.b2);
    b3_2 = double(sol2.b3);
    
    % For clarity, define a function handle for the diverging section:
    divergingFun = @(x) b0_2 + b1_2*x + b2_2*x.^2 + b3_2*x.^3;
end


