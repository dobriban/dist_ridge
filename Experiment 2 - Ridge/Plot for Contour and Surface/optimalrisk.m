function a  = optimalrisk(alpha,gamma)

a = gamma.*stieltjes_mp(-gamma./(alpha.^2),gamma);
