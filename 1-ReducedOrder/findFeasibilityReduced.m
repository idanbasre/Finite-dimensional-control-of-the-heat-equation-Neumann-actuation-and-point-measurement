function [topt, P_opt, alpha1_opt] ...
= findFeasibilityReduced(N0, N, q, delta, K0_hat, L, F0)

    lambda_NP1 = pi^2*N^2;
    gamma = 2;
    mu_NP1 = 1 + gamma + lambda_NP1/gamma;
    
    % Setting LMIs
    setlmis([]);
    P = lmivar(1,[2*N0+1,1]);
    alpha1 = lmivar(1, [1,1]);

    % Psi_2 < 0
    lmiterm([1 1 1 P], 1, F0, 's');
    lmiterm([1 1 1 P], 2*delta, 1);
    lmiterm([1 1 1 alpha1], 4*(K0_hat'*K0_hat)/(N*pi^2), 1);    

    lmiterm([1 1 2 P], 1, L);
    
    lmiterm([1 2 2 0], -2*(lambda_NP1^2-(q+delta)*lambda_NP1)/mu_NP1);
    
    lmiterm([1 2 3 0], 1);
    
    lmiterm([1 3 3 alpha1], -mu_NP1/lambda_NP1, 1);
    
    % alpha1 > 0
    lmiterm([-2 1 1 alpha1], 1, 1);
    
    % P > 0
    lmiterm([-3 1 1 P], 1, 1);
    
    % LMI solving
    lmis = getlmis;
    [topt, xopt] = feasp(lmis,[0, 300, -5, 30, 0]);
    
    % Getting Descision Vars
    P_opt = dec2mat(lmis, xopt, P);
    alpha1_opt = dec2mat(lmis, xopt, alpha1);
end