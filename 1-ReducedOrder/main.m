%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Finite dimensional control of the heat equation          %%%
%%%             Neumann actuation and point measurement              %%%
%%%                   Reduced Order - H1 Stability                   %%%
%%%                                                                  %%%
%%%                                                                  %%%
%%%                                                                  %%%
%%%                     Author: Idan Basre                           %%%
%%%                       November 2021                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Initial Parameters
% ------------------
% Env
logs_folder = "Logs";

% System
q = 0;
N0 = 0;
delta = 70;

% Gains
delta_L0 =  10;
delta_K0 = 10;
% ------------------

% Calculating System
[A0, B0, B_ns, C0, C_ns, C_norm_sq] = getSystem(N0, q);

% Calculating Gains
[K0, L0] = calcGains(N0, A0, B0, C0, delta_L0, delta_K0);

% Solving LMIs   
t_opt = 2;
N = N0 + 1;
dN = 1;
while t_opt > 0 && N < 20
    [K, L, F, F1, F2, F3] = getClosedLoopSystemReduced(N0, A0, B0, C0, K0, L0);
    
    [t_opt, P_opt, alpha_opt, S0_opt, R0_opt, S1_opt, R1_opt, S2_opt, R2_opt, G1_opt, G2_opt] ...
        = findFeasibilityDelayedReduced(N0, N, q, tau_m, eta_m, r, K, L, F, F1, F2, F3, delta, delta_1, C_norm_sq, C_ns);
    
    log = PrepareLog(N0, N, dN, t_opt, P_opt, alpha_opt, S0_opt, R0_opt, S1_opt, R1_opt, S2_opt, R2_opt, G1_opt, G2_opt);
    N
    N = N+dN;
end

% Dumping logs
logs_path = logs_folder + "\" + datestr(now,'yyyy_mm_dd_HH_MM_SS');
xlswrite(logs_path, log);
save(logs_path);
