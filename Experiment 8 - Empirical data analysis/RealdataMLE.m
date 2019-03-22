%Real data experiment, comparison of theoretical ARE and emipirical ARE 
%MLE
%load default_features_1059_tracks.txt
%data = csvread('Frogs_MFCCs 2.csv',1,1);
%data = default_features_1059_tracks;
load YearPredictionMSD.txt;
data = YearPredictionMSD;
[n,p] = size(data);
X = data(:,2:p);
Y = data(:,1);
p = p-1;

n_array = [500,1000,2000,5000];

%randomly sample from X and Y

for w = 1:length(n_array)
rng(2);
savefigs=1;    closefigs=1;

n_test = n_array(w);
sample = randperm(n);
sample = sample(1:n_test);
sample = sort(sample);

X_test = X(sample,:);
X_test = X_test-sum(X_test)/n_test;
X_test = whiten(X_test,0.00001);
Y_test = Y(sample,:);
Y_test = Y_test-sum(Y_test)/n_test;
Y_test = Y_test/sqrt(Y_test'*Y_test);
m = 100; %number of machines
gamma = p/n_test; %aspect ratio
Y_pred = X_test*inv(X_test'*X_test)*X_test'*Y_test;
SSres = (Y_pred-Y_test)'*(Y_pred-Y_test);
SStot = (Y_test-sum(Y_test)/n_test)'*(Y_test-sum(Y_test)/n_test);
Rsquare = 1-SSres/SStot;

%Estimation of signal and noise
[U,S,~] = svd(X_test);
S_sq = diag(S).^2/p;
S_sq = [S_sq; zeros(n_test-p,1)];
tY = U'*Y_test;
V = tY.^2;
f = ones(n_test,2);
f(:,1) = S_sq;
%e = f\V;
%alpha = sqrt(e(1));
%sigma = sqrt(e(2));
%alpha = alpha/sigma;

fun = @(x) sum(log(f*x'))+sum(V./(f*x'))
x0 = [0.0001,0.0001];
options = optimset('MaxFunEvals',10000,'MaxIter',10000);
x = fminsearch(fun, x0, options);
alpha = sqrt(x(1)/x(2));

%Empirical ARE
ARES = zeros(1,m);
ARES(1) = 1;
for i = 2:m
    
    Xcell = cell(1,i);
    beta = zeros(p,i);
    Xcell{1} = X_test(1:floor(n_test/i),:);
    beta(:,1) = (Xcell{1}'*Xcell{1}+p/(alpha^2)*eye(p))\Xcell{1}'*Y_test(1:floor(n_test/i));
    Q = cell(1,i);
    Q{1} = (Xcell{1}'*Xcell{1}+p/(alpha^2)*eye(p))\(Xcell{1}'*Xcell{1});
    
    for j = 2:i-1
        Xcell{j} = X_test((j-1)*floor(n_test/i)+1:j*floor(n_test/i),:);
        Q{j} = (Xcell{j}'*Xcell{j}+p/(alpha^2)*eye(p))\(Xcell{j}'*Xcell{j});
        beta(:,j) = (Xcell{j}'*Xcell{j}+p/(alpha^2)*eye(p))\Xcell{j}'*...
           Y_test((j-1)*floor(n_test/i)+1:j*floor(n_test/i));
    end
    Xcell{i} = X_test((i-1)*floor(n_test/i)+1:n_test,:);
    Q{i} = (Xcell{i}'*Xcell{i}+p/(alpha^2)*eye(p))\(Xcell{i}'*Xcell{i});
    beta(:,i) = (Xcell{i}'*Xcell{i}+p/(alpha^2)*eye(p))\Xcell{i}'*...
        Y_test((i-1)*floor(n_test/i)+1:n_test);
    
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
   Q1 = (X_test'*X_test+p/(alpha^2)*eye(p))\(X_test'*X_test);
   v1 = alpha^2/p*trace(Q1);
   A1 = alpha^2/p*trace(Q1*Q1);
   R1 = trace((X_test'*X_test+p/(alpha^2)*eye(p))\Q1);
   M1 = alpha^2-v1^2/(A1+R1);
   ARES(i) = M1/(alpha^2-v'*inv(A+R)*v) ;
end

k = 1:m;
ARET = gamma(1).*stieltjes_mp(-gamma(1)/(alpha)^2,gamma(1))./(alpha)^2.*(1-k+...
    (alpha)^2./(gamma(1).*stieltjes_mp(-k.*gamma(1)/(alpha)^2,k.*gamma(1))));
plot(k,ARET,'-.','LineWidth',3)
hold on 
plot(k,ARES,':','LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('RE');
xlim([min(k),max(k)]);
ylim([0,1.2]);
set(gca,'fontsize',20)
legend('Theoretical','Empirical');
str = sprintf( 'SNR=%.2f-gamma=%.2f-Rsquared=%.2f',alpha^2,gamma,Rsquare);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './Ridge-Realdata-n=%d-gamma=%.3f-SNR=%.2f.png',...
            n_test,gamma,alpha^2);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end

end
