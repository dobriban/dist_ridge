%Plots of phi for different alpha
cd('/Users/shengyue/Desktop')

alpha = 0.5;

%rng(2)
x = [0:0.1:10];
y = alpha^2./(2*x).*(-x/alpha^2+x-1+sqrt((-x/alpha^2+x-1).^2+4*x.^2/alpha^2));
plot(x,y,'LineWidth',3)
xlabel('x');
ylabel('\phi(x)');
set(gca,'fontsize',20);
xlim([min(x),max(x)]);
str = sprintf( 'alpha=%.2f',alpha);
title(str);

savefigs = 1;
  if savefigs==1
        filename = sprintf( './phi-plot-alpha=%.2f.png',alpha);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
  end
 