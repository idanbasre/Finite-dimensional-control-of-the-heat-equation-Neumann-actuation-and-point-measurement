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
q = 1;
N0 = 1;
delta = 2;

% Gains
delta_L0 =  delta;
delta_K0 = delta;
% ------------------

% Calculating System
[A0, A0_hat, B0, B_ns, B0_tilde, C0, C_ns] = getSystem(N0, q);

% Calculating Gains
% [K0, L0] = calcGains(N0, A0, A0_hat, B0_tilda, C0, delta_L0, delta_K0);
L0 = 3.5;
K0 = [1, 3.5];

% Solving LMIs   
t_opt = 2;
N = N0 + 1;
dN = 1;
while t_opt > 0 && N < 15
    [K0_hat, L0_tilde, F0, L] = getClosedLoopSystemReduced(N0, q, A0, A0_hat, B0_tilde, C0, K0, L0);
    [topt, P_opt, alpha1_opt] = findFeasibilityReduced(N0, N, q, delta, K0_hat, L, F0);
    log = prepareLog(N0, N, dN, q, delta, K0_hat, L, F0, A0_hat, A0, B0_tilde, C0, K0, L0, t_opt, P_opt, alpha1_opt);
    N
    N = N+dN;
end

% Dumping logs
logs_path = logs_folder + "\" + datestr(now,'yyyy_mm_dd_HH_MM_SS');
xlswrite(logs_path, log);
save(logs_path);


% T = cell2table(a(2:end,:),'VariableNames',a(1,:));
% writetable(T,logs_path);
