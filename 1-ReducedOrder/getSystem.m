function [A0, A0_hat, B0, B_ns, B0_tilde, C0, C_ns] = getSystem(N0, q)
    phi_n = @(x,n) sqrt(2)*cos(pi*(n-1)*x);
    N = 30;
    
    % A0
    A0 = diag(-pi^2*((1:N0)-1).^ 2 + q);
    
    %B0
    B_ns = zeros(N,1);
    B_ns(1,1) = 1/3;
    for n = 2: N
        B_ns(n,1) = -sqrt(2)*(n-1)^(-2)*pi^(-2);
    end
    B0 = B_ns(1:N0);
    
    assert(all(B_ns));
    
    %A0_hat
    A0_hat = [ 0, zeros(1, N0);
               q*B0,   A0     ];
    
    % B0_tilda
    B0_tilde = [1; -B0];

    %C0  
    C_n = @(n) phi_n(0,n);
    C_ns = zeros(1,N);
    C_ns(1,1) = 1; % phi_1(0)=1
    for n = 2: N
        C_ns(1,n) = C_n(n);
    end
    C0 = C_ns(1:N0);
    
    assert(all(C0));
end