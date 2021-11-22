function [K0, L0] = calcGains(N0, A0, B0, C0, delta_L0, delta_K0)

    % L0
    setlmis([]);
    P = lmivar(1, [N0+1,1]);
    Y = lmivar(2, [N0+1,1]); % Y = P*L0

    lmiterm([1 1 1 P], A0', 1, 's');
    lmiterm([1 1 1 Y], 1, -C0, 's');
    lmiterm([1 1 1 P], 2*delta_L0, 1);

    lmiterm([-2 1 1 P], 1, 1);

    lmis = getlmis;
    [topt, xopt] = feasp(lmis,[0, 0, -1, 0, 1]);
    assert(topt < 0)
    Y_opt = dec2mat(lmis, xopt, Y);
    P_opt = dec2mat(lmis, xopt, P);
    L0 = P_opt\Y_opt;

    % Tests L0
    assert(max(eig(A0-L0*C0)) < 0);

    % K0
    setlmis([]);
    Pinv = lmivar(1, [N0+1,1]);
    Q = lmivar(2, [1, N0+1]); % Q = K0*inv(P)

    lmiterm([1 1 1 Pinv], A0, 1, 's');
    lmiterm([1 1 1 Q], B0, 1, 's');
    lmiterm([1 1 1 Pinv], 2*delta_K0, 1); 

    lmiterm([-2 1 1 Pinv], 1, 1);

    lmiterm([3 1 1 Pinv], 1, 1);
    lmiterm([-3 1 1 0], 1.0e8);

    lmis = getlmis;
    [topt, xopt] = feasp(lmis,[0, 0, -1, 0, 1]);
    assert(topt < 0)
    Q_opt = dec2mat(lmis, xopt, Q);
    Pinv_opt = dec2mat(lmis, xopt, Pinv);
    K0 = Q_opt/Pinv_opt;

    % Tests K0
    assert(max(eig(A0+B0*K0)) < 0);
end
