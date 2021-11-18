function log = prepareLog(N0, N, dN, t_opt, P_opt, alpha_opt, S0_opt, R0_opt, S1_opt, R1_opt, S2_opt, R2_opt, G1_opt, G2_opt);
    if N == N0 + dN
        % Results Legend
        log(1,1) = "topt";
    end

    % Adding Current Iteration to res Matrix
    log(N,1) = t_opt;
end