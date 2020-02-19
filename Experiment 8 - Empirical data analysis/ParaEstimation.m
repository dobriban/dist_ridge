load YearPredictionMSD.txt;
data = YearPredictionMSD;
[n,p] = size(data);
Xtrain = data(1:100000,2:p);
Ytrain = data(1:100000,1);
Xtest = data(463716:n,2:p);
Ytest = data(463716:n,1);
p = p-1; %number of features

rng(0);

n_train = 10000;
n_test = 1000;

Trainsample = randperm(size(Xtrain,1));
Trainsample = Trainsample(1:n_train);
Trainsample = sort(Trainsample);
X_train = Xtrain(Trainsample,:);
Y_train = Ytrain(Trainsample,:);

Testsample = randperm(size(Xtest,1));
Testsample = Testsample(1:n_test);
Testsample = sort(Testsample);
X_test = Xtest(Testsample,:);
Y_test = Ytest(Testsample,:);

%Normalize data
X_train = X_train-mean(X_train);
X_train = normc(X_train);
Ymean = sum(Y_train)/n_train;
Ynorm = norm(Y_train);
Y_train = Y_train-Ymean;
Y_train = Y_train/Ynorm;
X_test = X_test-mean(X_test);
X_test = normc(X_test);

%Estimation of signal and noise by MLE
[U,S,~] = svd(X_train);
S_sq = [diag(S).^2/p; zeros(n_train-p,1)];
tY = U'*Y_train;
V = tY.^2;
f = ones(n_train,2);
f(:,1) = S_sq;
fun = @(x) sum(log(f*x'))+sum(V./(f*x'));
x0 = [0.0001,0];
options = optimset('MaxFunEvals',400,'MaxIter', 400);
x = fminsearch(fun, x0, options);
alpha = sqrt(x(1)/x(2));