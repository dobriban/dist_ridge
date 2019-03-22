function w  = opt_w(k,gamma,alpha)

lambda = k*gamma./alpha.^2;
m  = stieltjes_mp(-lambda,k*gamma);
num  = 1./(lambda.*m);
denom = 1-k+k*num;
w = num./denom;


