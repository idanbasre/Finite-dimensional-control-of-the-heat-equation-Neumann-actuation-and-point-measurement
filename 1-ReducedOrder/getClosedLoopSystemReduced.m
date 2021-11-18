function [K, L, F, F1, F2, F3] = getClosedLoopSystemReduced(N0, A0, B0, C0, K0, L0)
    K = [K0, zeros(1, N0+1)];
    L = [ L0;
         -L0];
    
    F = [     A0,                L0*C0  ;
          zeros(N0+1,N0+1),    A0-L0*C0  ];

    F1 = [   zeros(N0+1),        L0*C0    ;
             zeros(N0+1),       -L0*C0   ];

    F2 = [   B0;
          zeros(N0+1,1);];
      
    F3 = [zeros(N0+1,1)
             B0;];

end