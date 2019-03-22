cd('/Users/shengyue/Desktop')
alpha = [0.5,1,2,5];

savefigs=1;    closefigs=1;

x = 0:0.1:10;
y1 = x.*stieltjes_mp(-x/alpha(1)^2,x)./alpha(1)^2.*(1+...
    alpha(1)^2./(x*(1+alpha(1)^2)));
y1(1) = 1/(1+0.25);
plot(x,y1,'LineWidth',2)
hold on 
y2 = x.*stieltjes_mp(-x/alpha(2)^2,x)./alpha(2)^2.*(1+...
    alpha(2)^2./(x*(1+alpha(2)^2)));
y2(1) = 0.5;
plot(x,y2,'LineWidth',2)
y3 = x.*stieltjes_mp(-x/alpha(3)^2,x)./alpha(3)^2.*(1+...
    alpha(3)^2./(x*(1+alpha(3)^2)));
y3(1) = 1/(1+4);
plot(x,y3,'LineWidth',2)
y4 = x.*stieltjes_mp(-x/alpha(4)^2,x)./alpha(4)^2.*(1+...
    alpha(4)^2./(x*(1+alpha(4)^2)));
y4(1) = 1/(1+25);
plot(x,y4,'LineWidth',2)
hold off
xlabel('gamma');
ylabel('limit of ARE');
set(gca,'fontsize',20)
legend('alpha=0.5','alpha=1','alpha=2','alpha=5');

  if savefigs==1
        filename = ...
            sprintf( './limit-of-ARE.png');
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
  end
