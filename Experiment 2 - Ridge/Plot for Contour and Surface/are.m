function a  = are(alpha,gamma)

A1 = gamma.*stieltjes_mp(-gamma./(alpha.^2),gamma)./(alpha.^2);
A2=  1+alpha.^2./(gamma.*(1+alpha.^2));
a = A1.*A2;
