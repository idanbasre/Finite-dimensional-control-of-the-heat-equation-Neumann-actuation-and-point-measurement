function log = prepareLog(N0, N, dN, q, delta, K0_hat, L, F0, A0_hat, A0, B0_tilde, C0, K0, L0, t_opt, P_opt, alpha1_opt)
    if N == N0 + dN
        % Results Legend
        log(1,1) = "topt";
    end

    % Adding Current Iteration to res Matrix
    log(N,1) = t_opt;
end