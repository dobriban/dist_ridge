function m  = stieltjes_mp(lambda,gamma)

a = lambda+gamma-1;
m = (a+(a.^2-4.*gamma.*lambda).^(1/2))./(-2.*gamma.*lambda);
