function [A0, B0, B_ns, C0, C_ns, C_norm_sq] = getSystem(N0, q)
    phi_n = @(x,n) sqrt(2)*cos(pi*n*x);
    N = 30;
    
    % A0
    A0 = diag(-pi^2*((0:N0)).^ 2 + q);
    
    %B0
    B_ns = zeros(N,1);
    B_ns(1,1) = 1;
    for n = 2: N
        B_ns(n,1) = sqrt(2)*(-1)^(n-1);
    end
    B0 = B_ns(1:N0+1);
    
    assert(all(B_ns));

    %C0
    C = @(x) heaviside(x-0.3) - heaviside(x-0.9);
    C_n = @(n) integral(@(x) C(x).*phi_n(x,n),0,1);
    C_norm_sq = integral(@(x) C(x).*C(x),0,1);
    C_ns = zeros(1,N);
    C_ns(1,1) = 0.6;
    for n = 2: N
        C_ns(1,n) = C_n(n-1);
    end
    C0 = C_ns(1:N0+1);
    C_norm_N = C_norm_sq - norm(C_ns(1:6+1))^2;
    
    assert(all(C0));
end