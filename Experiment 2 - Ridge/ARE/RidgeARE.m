cd('/Users/shengyue/Desktop')
n = 1000;
p = [200,500,1000,2000,5000];
alpha_arr =[0.5,1,2,5];
gamma = p/n;
L = length(alpha_arr);

for l=1:L
    
alpha = alpha_arr(l);

savefigs=1;    closefigs=1;

k = 1:0.1:20;
ARE1 = gamma(1).*stieltjes_mp(-gamma(1)/alpha^2,gamma(1))./alpha^2.*(1-k+...
    alpha^2./(gamma(1).*stieltjes_mp(-k.*gamma(1)/alpha^2,k.*gamma(1))));
plot(k,ARE1,'LineWidth',3)
hold on 
ARE2 = gamma(2).*stieltjes_mp(-gamma(2)/alpha^2,gamma(2))./alpha^2.*(1-k+...
    alpha^2./(gamma(2).*stieltjes_mp(-k.*gamma(2)/alpha^2,k.*gamma(2))));
plot(k,ARE2,'LineWidth',3)
ARE3 = gamma(3).*stieltjes_mp(-gamma(3)/alpha^2,gamma(3))./alpha^2.*(1-k+...
    alpha^2./(gamma(3).*stieltjes_mp(-k.*gamma(3)/alpha^2,k.*gamma(3))));
plot(k,ARE3,'LineWidth',3)
ARE4 = gamma(4).*stieltjes_mp(-gamma(4)/alpha^2,gamma(4))./alpha^2.*(1-k+...
    alpha^2./(gamma(4).*stieltjes_mp(-k.*gamma(4)/alpha^2,k.*gamma(4))));
plot(k,ARE4,'LineWidth',3)
ARE5 = gamma(5).*stieltjes_mp(-gamma(5)/alpha^2,gamma(5))./alpha^2.*(1-k+...
    alpha^2./(gamma(5).*stieltjes_mp(-k.*gamma(5)/alpha^2,k.*gamma(5))));
plot(k,ARE5,'LineWidth',3)
hold off
xlabel('Number of Machines');
ylabel('ARE');
xlim([min(k),max(k)]);
set(gca,'fontsize',20)
legend('gamma=0.2','gamma=0.5','gamma=1','gamma=2','gamma=5');
str = sprintf( 'alpha=%.2f',alpha);
title(str);

  if savefigs==1
        filename = ...
            sprintf( './ARE-alpha=%.2f.png',alpha);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
end
