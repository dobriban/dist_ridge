function a  = oe(alpha,gamma)

A1 = gamma.*stieltjes_mp(-gamma./(alpha.^2),gamma);
A2 = alpha.^2./(1+alpha.^2./(gamma.*(1+alpha.^2)));
a = (1+A1)./(1+A2);