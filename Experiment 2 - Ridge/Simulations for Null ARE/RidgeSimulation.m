%Simulations of ARE for ridge regression, null case
cd('/Users/shengyue/Desktop')
n = 1000;
p_arr = [100,500];
alpha = 1;
L = length(p_arr);
m = 20;

for l=1:L
    
p = p_arr(l);
gamma = p/n;
savefigs=1;    closefigs=1;
rng(3);
X = normrnd(0,1,[n,p]);
beta = normrnd(0,alpha/sqrt(p),[p,1]);

ARES = zeros(1,m);
ARES(1)=1;
for i = 2:m
    Xcell = cell(1,i);
    Xcell{1} = X(1:floor(n/i),:);
    Q = cell(1,i);
    Q{1} = (Xcell{1}'*Xcell{1}+p/(alpha^2)*eye(p))\(Xcell{1}'*Xcell{1});
    for j = 2:i-1
        Xcell{j} = X((j-1)*floor(n/i)+1:j*floor(n/i),:);
        Q{j} = (Xcell{j}'*Xcell{j}+p/(alpha^2)*eye(p))\(Xcell{j}'*Xcell{j});
    end
    Xcell{i} = X((i-1)*floor(n/i)+1:n,:);
    Q{i} = (Xcell{i}'*Xcell{i}+p/(alpha^2)*eye(p))\(Xcell{i}'*Xcell{i});
    
    B = zeros(p,i);
    for t = 1:i
        B(1:p,t) = Q{t}*beta;
    end
    v = zeros(i,1);
    for t = 1:i
        v(t) = trace((Xcell{t}'*Xcell{t}+p/(alpha^2)*eye(p))\Q{t});
    end
    A = diag(v);
    
M =beta'*(eye(p)-B*inv(B'*B+A)*(B'))*beta;
ARES(i) = (p^2/alpha^4*beta'*inv((X'*X+p/alpha^2*eye(p))^2)*beta+...
    trace(inv((X'*X+p/alpha^2*eye(p))^2)*(X'*X)))/M;
end
    


k = 1:m;
ARET = gamma(1).*stieltjes_mp(-gamma(1)/alpha^2,gamma(1))./alpha^2.*(1-k+...
    alpha^2./(gamma(1).*stieltjes_mp(-k.*gamma(1)/alpha^2,k.*gamma(1))));
plot(k,ARET,'LineWidth',3)
hold on 
plot(k,ARES,':','LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('ARE');
xlim([min(k),max(k)]);
ylim([0,1]);
set(gca,'fontsize',20)
legend('Location','southwest','Theoretical','Numerical');
str = sprintf( 'alpha=%.2f-gamma=%.2f-m=%.d',alpha,gamma,m);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './ARE-alpha=%.2f-gamma=%.2f-m=%.d.png',alpha,gamma,m);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end
