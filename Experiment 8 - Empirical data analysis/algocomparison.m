%Real data experiment, comparison of different algorithms

%load YearPredictionMSD.txt;
data = YearPredictionMSD;
[n,p] = size(data);
Xtrain = data(1:463715,2:p);
Ytrain = data(1:463715,1);
Xtest = data(463716:n,2:p);
Ytest = data(463716:n,1);
p = p-1;

m_array = [10,20,50,100,500];

T = 10;

MSE1 = zeros(length(m_array),T);
MSE2 = zeros(length(m_array),T);
MSE3 = zeros(length(m_array),T);

rng(0);

for q = 1:T

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

%Whiten data
X_train = X_train-sum(X_train)/n_train;
X_train = whiten(X_train,0.00001);
Ymean = sum(Y_train)/n_train;
Y_train = Y_train-Ymean;

X_test = X_test-sum(X_test)/n_test;
X_test = whiten(X_test,0.00001);
%Y_test = Y_test-sum(Y_test)/n_test;

%Estimation of signal and noise
[U,S,~] = svd(X_train);
S_sq = diag(S).^2/p;
S_sq = [S_sq; zeros(n_train-p,1)];
tY = U'*Y_train;
V = tY.^2;
f = ones(n_train,2);
f(:,1) = S_sq;
%e = f\V;
%alpha = sqrt(e(1));
%sigma = sqrt(e(2));
%alpha = alpha/sigma;

fun = @(x) sum(log(f*x'))+sum(V./(f*x'))
x0 = [0.0001,0.0001];
options = optimset('MaxFunEvals',1000,'MaxIter',1000);
x = fminsearch(fun, x0, options);
alpha = sqrt(x(1)/x(2));


for k = 1:length(m_array)
    
    i = m_array(k);
    Xcell = cell(1,i);
    beta = zeros(p,i);
    Beta = zeros(p,i);
    Xcell{1} = X_train(1:(n_train/i),:);
    beta(:,1) = (Xcell{1}'*Xcell{1}+p/(alpha^2)*eye(p))\Xcell{1}'*Y_train(1:(n_train/i));
    Beta(:,1) = (Xcell{1}'*Xcell{1}+p/(i*alpha^2)*eye(p))\Xcell{1}'*Y_train(1:(n_train/i));
    Q = cell(1,i);
    Q{1} = (Xcell{1}'*Xcell{1}+p/(alpha^2)*eye(p))\(Xcell{1}'*Xcell{1});
    
    for j = 2:i-1
        Xcell{j} = X_train((j-1)*(n_train/i)+1:j*(n_train/i),:);
        Q{j} = (Xcell{j}'*Xcell{j}+p/(alpha^2)*eye(p))\(Xcell{j}'*Xcell{j});
        beta(:,j) = (Xcell{j}'*Xcell{j}+p/(alpha^2)*eye(p))\Xcell{j}'*...
           Y_train((j-1)*(n_train/i)+1:j*(n_train/i));
        Beta(:,j) = (Xcell{j}'*Xcell{j}+p/(i*alpha^2)*eye(p))\Xcell{j}'*...
           Y_train((j-1)*(n_train/i)+1:j*(n_train/i));
    end
    
    Xcell{i} = X_train((i-1)*(n_train/i)+1:n_train,:);
    Q{i} = (Xcell{i}'*Xcell{i}+p/(alpha^2)*eye(p))\(Xcell{i}'*Xcell{i});
    beta(:,i) = (Xcell{i}'*Xcell{i}+p/(alpha^2)*eye(p))\Xcell{i}'*...
        Y_train((i-1)*(n_train/i)+1:n_train);
    Beta(:,i) = (Xcell{i}'*Xcell{i}+p/(i*alpha^2)*eye(p))\Xcell{i}'*...
        Y_train((i-1)*(n_train/i)+1:n_train);
    
    %Estimation of optimal weights
    v = zeros(i,1);
    A = zeros(i,i);
    R = zeros(i,1);
    for t = 1:i
        v(t) = alpha^2/p*trace(Q{t});
        R(t) = trace((Xcell{t}'*Xcell{t}+p/(alpha^2)*eye(p))\Q{t});
    end
    R = diag(R);
    
    for z = 1:i
        for w = 1:i
        A(w,z) = alpha^2/p*trace(Q{w}*Q{z});
        end
    end

%beta1 = beta*((A+R)\v);
weight = alpha^2/(alpha^2*i+(1-i)*i*p/n_train*stieltjes_mp(-i*p/(n_train*(alpha)^2),i*p/n_train));
beta1 = beta*ones(i,1)*weight;
Y_predict = X_test*beta1+Ymean;
MSE1(k,q) = sum((Y_test-Y_predict).^2)/n_test;

beta2 = Beta*ones(i,1)/i;
Y_predict2 = X_test*beta2+Ymean;
MSE2(k,q) = sum((Y_test-Y_predict2).^2)/n_test;

beta3 = beta(:,1);
Y_predict3 = X_test*beta3+Ymean;
MSE3(k,q) = sum((Y_test-Y_predict3).^2)/n_test;
end

end

MSE1mean = sum(MSE1,2)/T;
error1 = sqrt(sum((MSE1-MSE1mean).^2,2)/T);
MSE2mean = sum(MSE2,2)/T;
error2 = sqrt(sum((MSE2-MSE2mean).^2,2)/T);
MSE3mean = sum(MSE3,2)/T;
error3 = sqrt(sum((MSE3-MSE3mean).^2,2)/T);

errorbar(m_array,MSE1mean,error1,'-.o','LineWidth',3,'MarkerSize',...
    10,'MarkerEdgeColor','red','MarkerFaceColor','white')
hold on 
errorbar(m_array,MSE2mean,error2,'--s','LineWidth',3,'MarkerSize',...
    10,'MarkerEdgeColor','green','MarkerFaceColor','white')
errorbar(m_array,MSE3mean,error3,':*','LineWidth',3,'MarkerSize',...
    10,'MarkerEdgeColor','black','MarkerFaceColor','white')
hold off
xlabel('Number of Machines');
ylabel('MSE of Prediction');
set(gca,'fontsize',20)
legend('Optimal average','Naive average','1/k data','Location','northwest');



