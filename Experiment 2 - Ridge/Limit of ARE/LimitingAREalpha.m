cd('/Users/shengyue/Desktop')
gamma = [0.1,1,1.1,2];

savefigs=1;    closefigs=1;

x = 0:0.1:100;
y1 = gamma(1).*stieltjes_mp(-gamma(1)./(x.^2),gamma(1))./(x.^2).*(1+...
    (x.^2)./(gamma(1).*(1+x.^2)));
y1(1)=1;
plot(x,y1,'LineWidth',2)
hold on 
y2 = gamma(2).*stieltjes_mp(-gamma(2)./(x.^2),gamma(2))./(x.^2).*(1+...
    (x.^2)./(gamma(2).*(1+x.^2)));
y2(1)=1;
plot(x,y2,'LineWidth',2)
y3 = gamma(3).*stieltjes_mp(-gamma(3)./(x.^2),gamma(3))./(x.^2).*(1+...
    (x.^2)./(gamma(3).*(1+x.^2)));
y3(1)=1;
plot(x,y3,'LineWidth',2)
y4 = gamma(4).*stieltjes_mp(-gamma(4)./(x.^2),gamma(4))./(x.^2).*(1+...
    (x.^2)./(gamma(4).*(1+x.^2)));
y4(1)=1;
plot(x,y4,'LineWidth',2)
hold off
xlabel('alpha');
ylabel('limit of ARE');
set(gca,'fontsize',20)
legend('gamma=0.1','gamma=1','gamma=1.1','gamma=2');

  if savefigs==1
        filename = ...
            sprintf( './limit-of-ARE-alpha.png');
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
