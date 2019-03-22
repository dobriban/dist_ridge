%Simulations  with ridge regression
cd('C:\Dropbox\Projects\Distributed Ridge\Experiments\Experiment 4 - Ridge Regularization Param')
addpath('C:\Dropbox\Projects\Distributed Ridge\\Experiments\Code')
n = 3000;
p = 500;
alpha = 1;
gamma = p/n;
%%
for k=1:10
    
    %% generate the data
    T = 20;
    MSE_overall = zeros(T,1);
    MSE_naive = zeros(T,1);
    
    rng(2);
    n_mc = 50;
    for ii = 1:n_mc
        
        X = normrnd(0,1,[n,p]);
        beta = normrnd(0,alpha/sqrt(p),[p,1]);
        ep = normrnd(0,1,[n,1]);
        Y = X*beta+ ep;
        
        Xcell = cell(1,k);
        Xcell{1} = X(1:floor(n/k),:);
        Ycell = cell(1,k);
        Ycell{1} = Y(1:floor(n/k));
        for j = 2:k-1
            Xcell{j} = X((j-1)*floor(n/k)+1:j*floor(n/k),:);
            Ycell{j} = Y((j-1)*floor(n/k)+1:j*floor(n/k));
        end
        Xcell{k} = X((k-1)*floor(n/k)+1:n,:);
        Ycell{k} = Y((k-1)*floor(n/k)+1:n);
        
        %% search over lambdas
        T = 20;
        lambda_fac = linspace(0.1,2,T);
        ni_lambda_grid = p/(alpha^2)*lambda_fac;
        
        MSE_optw = zeros(T,1);
        MSE = zeros(T,1);
        for i=1:T
            ni_lambda = ni_lambda_grid(i);
            B = zeros(p,k);
            for t = 1:k
                B(:,t) = (Xcell{t}'*Xcell{t}+ni_lambda*eye(p))\(Xcell{t}'*Ycell{t});
            end
            hbeta = mean(B,2);
            MSE(i) = norm(hbeta-beta)^2;
            
            w  = opt_w(k,gamma,alpha);
            hbeta = w*k*mean(B,2);
            MSE_optw(i) = norm(hbeta-beta)^2;
        end
        
        MSE_naive = MSE_naive + MSE;
        MSE_overall = MSE_overall + MSE_optw;
    end
    MSE_mean_n = MSE_naive/n_mc;
    MSE_mean = MSE_overall/n_mc;
    
    %%
    rng(2);
    savefigs=1;    closefigs=0; a = {'-','--','-.',':'};
    figure, hold on
    h1 = plot(lambda_fac,MSE_mean,'LineWidth',3,'color',rand(1,3));
    set(h1,'LineStyle',a{1});
    pa=1; %your point goes here
    x=[pa,pa];
    y=get(gca,'Ylim');
    h2 = plot(x,y,'linewidth',3,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
    h3 = plot(lambda_fac,MSE_mean_n,'LineWidth',3,'color',rand(1,3));
    set(h3,'LineStyle',a{3});
    
    xlabel('lambda');
    ylabel('MSE');
    xlim([min(lambda_fac),max(lambda_fac)]);
    ylim([min(MSE_mean),max(MSE_mean_n)]);
    set(gca,'fontsize',20)
    legend([h1,h2,h3],{'MSE opt','theo opt lambda', 'MSE avg'},'location','Best')
    %legend('Location','southwest','Theoretical','Numerical');
    %str = sprintf( 'alpha=%.2f-gamma=%.2f-m=%.d',alpha,gamma,m);
    %title(str);
    
    if savefigs==1
        filename = ...
            sprintf( './MSE-fn-lambda-alpha=%.2f-gamma=%.2f-k=%.d-n_mc=%d.png',alpha,gamma,k,n_mc);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end
end
