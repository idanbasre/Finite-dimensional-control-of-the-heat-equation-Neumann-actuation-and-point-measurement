function [K0_hat, L0_tilde, F0, L] = getClosedLoopSystemReduced(N0, q, A0, A0_hat, B0_tilde, C0, K0, L0)
    L0_tilde = [0; L0];
    K0_hat = [K0-[q, zeros(N0,1)], zeros(N0,1)];
    
    F0 = [   A0_hat+B0_tilde*K0,    L0_tilde*C0  ;
              zeros(N0,N0+1),         A0-L0*C0  ];
          
    L = [ L0_tilde; -L0];

end