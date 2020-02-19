%Real data experiment, comparison of different algorithms, general cov

load YearPredictionMSD.txt;
data = YearPredictionMSD;
[n,p] = size(data);
Xtrain = data(1:463715,2:p);
Ytrain = data(1:463715,1);
Xtest = data(463716:n,2:p);
Ytest = data(463716:n,1);
p = p-1; %number of features
m_array = [1,2,5,10,20,50,100,200,500,1000]; %number of machines
T = 10; %number of trials

MSE1 = zeros(length(m_array),T); %MSE for optimal weighted average
MSE2 = zeros(length(m_array),T); %MSE for naive average
MSE3 = zeros(length(m_array),T); %MSE for 1/k data

rng(1);

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
%[U,S,~] = svd(X_train);
%S_sq = [diag(S).^2/p; zeros(n_train-p,1)];
%tY = U'*Y_train;
%V = tY.^2;
%f = ones(n_train,2);
%f(:,1) = S_sq;
%fun = @(x) sum(log(f*x'))+sum(V./(f*x'));
%x0 = [0.0001,0];
%options = optimset('MaxFunEvals',400,'MaxIter', 400);
%x = fminsearch(fun, x0, options);
%alpha = sqrt(x(1)/x(2));
alpha = 1.2247;

for k = 1:length(m_array)
    
    i = m_array(k);
    Xcell = cell(1,i);
    beta = zeros(p,i);
    optimal_weight = zeros(i,1);
    Beta = zeros(p,i);
    Xcell{1} = X_train(1:(n_train/i),:);
    beta(:,1) = (Xcell{1}'*Xcell{1}+p/alpha^2*eye(p))\Xcell{1}'*Y_train(1:(n_train/i));
    Beta(:,1) = (Xcell{1}'*Xcell{1}+p/(i*alpha^2)*eye(p))\Xcell{1}'*Y_train(1:(n_train/i));
    lambda = 3*i*p/(n_train*alpha^2);
    m = trace(inv(Xcell{1}'*Xcell{1}*i/n_train+lambda*eye(p)))/p;
    m_prime = trace((inv(Xcell{1}'*Xcell{1}*i/n_train+lambda*eye(p)))^2)/p;
    f_func = alpha^2*i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime)+i*p/n_train*(m-lambda*m_prime);
    g_func = alpha^2*(1-2*lambda*m+lambda^2*m_prime-i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime));
    optimal_weight(1) = alpha^2*(1-lambda*m)/(f_func+i*g_func);
    
    for j = 2:i-1
        Xcell{j} = X_train((j-1)*(n_train/i)+1:j*(n_train/i),:);
        beta(:,j) = (Xcell{j}'*Xcell{j}+p/alpha^2*eye(p))\Xcell{j}'*...
           Y_train((j-1)*(n_train/i)+1:j*(n_train/i));
        Beta(:,j) = (Xcell{j}'*Xcell{j}+p/(i*alpha^2)*eye(p))\Xcell{j}'*...
           Y_train((j-1)*(n_train/i)+1:j*(n_train/i));
       
        m = trace(inv(Xcell{j}'*Xcell{j}*i/n_train+lambda*eye(p)))/p;
        m_prime = trace((inv(Xcell{j}'*Xcell{j}*i/n_train+lambda*eye(p)))^2)/p;
        f_func = alpha^2*i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime)+i*p/n_train*(m-lambda*m_prime);
        g_func = alpha^2*(1-2*lambda*m+lambda^2*m_prime-i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime));
        optimal_weight(j) = alpha^2*(1-lambda*m)/(f_func+i*g_func);
    end
    
    Xcell{i} = X_train((i-1)*(n_train/i)+1:n_train,:);
    beta(:,i) = (Xcell{i}'*Xcell{i}+p/alpha^2*eye(p))\Xcell{i}'*...
        Y_train((i-1)*(n_train/i)+1:n_train);
    Beta(:,i) = (Xcell{i}'*Xcell{i}+p/(i*alpha^2)*eye(p))\Xcell{i}'*...
        Y_train((i-1)*(n_train/i)+1:n_train);
    
    m = trace(inv(Xcell{i}'*Xcell{i}*i/n_train+lambda*eye(p)))/p;
    m_prime = trace((inv(Xcell{i}'*Xcell{i}*i/n_train+lambda*eye(p)))^2)/p;
    f_func = alpha^2*i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime)+i*p/n_train*(m-lambda*m_prime);
    g_func = alpha^2*(1-2*lambda*m+lambda^2*m_prime-i*p/n_train*lambda^2*(m-lambda*m_prime)^2/(1-i*p/n_train+i*p/n_train*lambda*m_prime));
    optimal_weight(i) = alpha^2*(1-lambda*m)/(f_func+i*g_func);
    
beta1 = beta*optimal_weight;
Y_predict = (X_test*beta1)*Ynorm+Ymean;
MSE1(k,q) = sum((Y_test-Y_predict).^2)/n_test;

beta2 = Beta*ones(i,1)/i;
Y_predict2 = (X_test*beta2)*Ynorm+Ymean;
MSE2(k,q) = sum((Y_test-Y_predict2).^2)/n_test;

beta3 = beta(:,1);
Y_predict3 = (X_test*beta3)*Ynorm+Ymean;
MSE3(k,q) = sum((Y_test-Y_predict3).^2)/n_test;
end

end

MSE1mean = sum(MSE1,2)/T;
error1 = sqrt(sum((MSE1-MSE1mean).^2,2)/T)/5;
MSE2mean = sum(MSE2,2)/T;
error2 = sqrt(sum((MSE2-MSE2mean).^2,2)/T)/5;
MSE3mean = sum(MSE3,2)/T;
error3 = sqrt(sum((MSE3-MSE3mean).^2,2)/T)/5;

errorbar(m_array,MSE1mean,error1,'-.o','LineWidth',3,'MarkerSize',...
    5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
hold on 
errorbar(m_array,MSE2mean,error2,'--s','LineWidth',3,'MarkerSize',...
    5,'MarkerEdgeColor','red','MarkerFaceColor','red')
errorbar(m_array,MSE3mean,error3,':*','LineWidth',3,'MarkerSize',...
    5,'MarkerEdgeColor','yellow','MarkerFaceColor','yellow')
hold off
xlabel('Number of Machines');
ylabel('MSE of Prediction');
ylim([105 122]);
set(gca,'fontsize',20)
legend('Optimal weighted average','Naive average','1/k data','Location','southeast');