function [topt, P_opt, alpha_opt, S0_opt, R0_opt, S1_opt, R1_opt, S2_opt, R2_opt, G1_opt, G2_opt] ...
= findFeasibilityReduced(N0, N, q, tau_m, eta_m, r, K, L, F, F1, F2, F3, delta_0, delta_1, C_norm_sq, C_ns)

    delta = delta_0 - delta_1;
    C_norm_N_sq = C_norm_sq - norm(C_ns(1:N+1))^2;
    lambda_NP1 = pi^2*(N+1)^2;
    
    % Setting LMIs
    setlmis([]);
    P = lmivar(1,[2*(N0+1),1]);
    S0 = lmivar(1, [1,1]);
    R0 = lmivar(1, [1,1]);
    S1 = lmivar(1, [1,1]);
    R1 = lmivar(1, [1,1]);
    G1 = lmivar(1, [1,1]);
    S2 = lmivar(1,[2*(N0+1),1]);
    R2 = lmivar(1,[2*(N0+1),1]);
    G2 = lmivar(1,[2*(N0+1),1]);
    alpha = lmivar(1, [1,1]);

    % Psi_2 < 0
    % Omega_1
    lmiterm([1 1 1 P], 1, F, 's');
    lmiterm([1 1 1 P], 2*delta, 1);
    lmiterm([1 1 1 R0], exp(-2*delta_0*r)*K', K);
    lmiterm([1 1 1 S2], 1-exp(-2*delta_0*tau_m), 1);
    lmiterm([1 1 1 S0], K', K);    

    lmiterm([1 1 2 P], 1, F2*K);
    lmiterm([1 1 2 R0], -exp(-2*delta_0*r)*K', K);
    
    lmiterm([1 1 3 P], 1, L);
    
    lmiterm([1 2 2 alpha], 2/((pi^2)*(N)), K'*K);
    lmiterm([1 2 2 R0], exp(-2*delta_0*r)*K', K);
    lmiterm([1 2 2 S1], (exp(-2*delta_0*r)-exp(-2*delta_0*(r+eta_m)))*K', K);
    lmiterm([1 2 2 S0], -exp(-2*delta_0*r)*K', K);
    
    lmiterm([1 3 3 0], -2*delta_1/C_norm_N_sq);

    
    % Theta_1
    lmiterm([1 1 4 P], 2*delta_1, 1);
    lmiterm([1 1 4 P], -1, F1);
    lmiterm([1 1 4 S2], exp(-2*delta_0*tau_m), 1);
    
    lmiterm([1 1 5 S2], exp(-2*delta_0*tau_m), 1);
    
    
    % Theta_2
    lmiterm([1 1 6 P], -1, F3);
    lmiterm([1 1 6 S1], exp(-2*delta_0*(r+eta_m))*K', 1);
    
    lmiterm([1 1 7 S1], exp(-2*delta_0*(r+eta_m))*K', 1);
    
    
    % Omega_2
    lmiterm([1 4 4 P], -2*delta_1, 1);
    lmiterm([1 4 4 R2], -exp(-2*delta_0*tau_m), 1);
    lmiterm([1 4 4 S2], -exp(-2*delta_0*tau_m), 1);
   
    lmiterm([1 4 5 S2], -exp(-2*delta_0*tau_m), 1);
    lmiterm([1 4 5 G2], -exp(-2*delta_0*tau_m), 1);
    
    lmiterm([1 5 5 R2], -exp(-2*delta_0*tau_m), 1);
    lmiterm([1 5 5 S2], -exp(-2*delta_0*tau_m), 1);
    
    
    % Omega_3
    lmiterm([1 6 6 alpha], 2/((pi^2)*(N)), 1);
    lmiterm([1 6 6 R1], -exp(-2*delta_0*(r+eta_m)), 1);
    lmiterm([1 6 6 S1], -exp(-2*delta_0*(r+eta_m)), 1);
    
    lmiterm([1 6 7 S1], -exp(-2*delta_0*(r+eta_m)), 1);
    lmiterm([1 6 7 G1], -exp(-2*delta_0*(r+eta_m)), 1);
    
    lmiterm([1 7 7 R1], -exp(-2*delta_0*(r+eta_m)), 1);
    lmiterm([1 7 7 S1], -exp(-2*delta_0*(r+eta_m)), 1);
    
    
    % Lambda part
    lmiterm([1 1 1 R2], tau_m^2*F', F);
    lmiterm([1 1 1 R1], eta_m^2*F'*K', K*F);
    lmiterm([1 1 1 R0], r^2*F'*K', K*F);  
    
    lmiterm([1 1 2 R2], tau_m^2*F', F2*K);
    lmiterm([1 1 2 R1], eta_m^2*F'*K', K*F2*K);
    lmiterm([1 1 2 R0], r^2*F'*K', K*F2*K);
    
    lmiterm([1 1 3 R2], tau_m^2*F', L);
    lmiterm([1 1 3 R1], eta_m^2*F'*K', K*L);
    lmiterm([1 1 3 R0], r^2*F'*K', K*L);
    
    lmiterm([1 1 4 R2], -tau_m^2*F', F1);
    lmiterm([1 1 4 R1], -tau_m^2*F'*K', K*F1);
    lmiterm([1 1 4 R0], -r^2*F'*K', K*F1);
    
    lmiterm([1 1 6 R2], -tau_m^2*F', F3);
    lmiterm([1 1 6 R1], -tau_m^2*F'*K', K*F3);
    lmiterm([1 1 6 R0], -r^2*F'*K', K*F3);
    
    lmiterm([1 2 2 R2], tau_m^2*K'*F2', F2*K);
    lmiterm([1 2 2 R1], eta_m^2*K'*F2'*K', K*F2*K);
    lmiterm([1 2 2 R0], r^2*K'*F2'*K', K*F2*K);  
    
    lmiterm([1 2 3 R2], tau_m^2*K'*F2', L);
    lmiterm([1 2 3 R1], eta_m^2*K'*F2'*K', K*L);
    lmiterm([1 2 3 R0], r^2*K'*F2'*K', K*L);
    
    lmiterm([1 2 4 R2], -tau_m^2*K'*F2', F1);
    lmiterm([1 2 4 R1], -eta_m^2*K'*F2'*K', K*F1);
    lmiterm([1 2 4 R0], -r^2*K'*F2'*K', K*F1);
    
    lmiterm([1 2 6 R2], -tau_m^2*K'*F2', F3);
    lmiterm([1 2 6 R1], -eta_m^2*K'*F2'*K', K*F3);
    lmiterm([1 2 6 R0], -r^2*K'*F2'*K', K*F3);
    
    lmiterm([1 3 3 R2], tau_m^2*L', L);
    lmiterm([1 3 3 R1], eta_m^2*L'*K', K*L);
    lmiterm([1 3 3 R0], r^2*L'*K', K*L);
    
    lmiterm([1 3 4 R2], -tau_m^2*L', F1);
    lmiterm([1 3 4 R1], -eta_m^2*L'*K', K*F1);
    lmiterm([1 3 4 R0], -r^2*L'*K', K*F1);
    
    lmiterm([1 3 6 R2], -tau_m^2*L', F3);
    lmiterm([1 3 6 R1], -eta_m^2*L'*K', K*F3);
    lmiterm([1 3 6 R0], -r^2*L'*K', K*F3);
    
    lmiterm([1 4 4 R2], tau_m^2*F1', F1);
    lmiterm([1 4 4 R1], eta_m^2*F1'*K', K*F1);
    lmiterm([1 4 4 R0], r^2*F1'*K', K*F1);
    
    lmiterm([1 4 6 R2], tau_m^2*F1', F3);
    lmiterm([1 4 6 R1], eta_m^2*F1'*K', K*F3);
    lmiterm([1 4 6 R0], r^2*F1'*K', K*F3);
    
    lmiterm([1 6 6 R2], tau_m^2*F3', F3);
    lmiterm([1 6 6 R1], eta_m^2*F3'*K', K*F3);
    lmiterm([1 6 6 R0], r^2*F3'*K', K*F3);

    % W_n_2 < 0
    lmiterm([2 1 1 0], -lambda_NP1 + q + delta_0);

    lmiterm([2 1 2 0], sqrt(0.5));

    lmiterm([2 2 2 alpha], -1/lambda_NP1, 1);

    % P > 0
    lmiterm([-3 1 1 P], 1, 1);

    % S0 > 0
    lmiterm([-4 1 1 S0], 1, 1);
    
    % R0 > 0
    lmiterm([-5 1 1 R0], 1, 1);
    
    % S1 > 0
    lmiterm([-6 1 1 S1], 1, 1);
    
    % R1 > 0
    lmiterm([-7 1 1 R1], 1, 1);
    
    % S2 > 0
    lmiterm([-8 1 1 S2], 1, 1);
    
    % R2 > 0
    lmiterm([-9 1 1 R2], 1, 1);
    
    % alpha > 0
    lmiterm([-10 1 1 alpha], 1, 1);
    
    % LMI solving
    lmis = getlmis;
    [topt, xopt] = feasp(lmis,[0, 300, -5, 30, 0]);
    
    % Getting Descision Vars
    P_opt = dec2mat(lmis, xopt, P);
    alpha_opt = dec2mat(lmis, xopt, alpha);
    S0_opt = dec2mat(lmis, xopt, S0);
    R0_opt = dec2mat(lmis, xopt, R0);
    S1_opt = dec2mat(lmis, xopt, S1);
    R1_opt = dec2mat(lmis, xopt, R1);
    S2_opt = dec2mat(lmis, xopt, S2);
    R2_opt = dec2mat(lmis, xopt, R2);
    G1_opt = dec2mat(lmis, xopt, G1);
    G2_opt = dec2mat(lmis, xopt, G2);
end