function m  = stieltjes_prime(lambda,gamma)

a = lambda+gamma-1;
m = ((gamma-1).*(a.^2-4.*gamma.*lambda).^(1/2)+(gamma-1).^2-...
    (1+gamma).*lambda)./(2.*gamma.*(lambda).^2.*(a.^2-4.*gamma.*lambda).^(1/2));