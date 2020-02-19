%Simulations  with ridge regression
%cd('C:\Dropbox\Projects\Distributed Ridge\Experiments\Experiment 7 - Lambda Gen Cov')
n = 3000;
p = 500;
alpha = 1;
gamma = p/n;
design = 0;
switch design
    case 0
        %Design 1: Id
        rho = 0;
        top = rho.^(0:1:p-1);
        Sigma = toeplitz(top);
        Sigma= diag(eig(Sigma));
        T_r = Sigma^(1/2);
        param = rho;
    case 1
        %Design 1: Toeplitz
        rho = 0.9;
        top = rho.^(0:1:p-1);
        Sigma = toeplitz(top);
        Sigma= diag(eig(Sigma));
        T_r = Sigma^(1/2);
        param = rho;
end
%%
k_arr = [1,2,5,10];
for ind=1:length(k_arr)
    k = k_arr(ind);
    
    T = 20;
    MSE_usual = zeros(T,1);
    
    rng(2);
    n_mc = 20;
    X = normrnd(0,1,[n,p])*T_r;
    Xcell = cell(1,k);
    Xcell{1} = X(1:floor(n/k),:);
    for j = 2:k-1
        Xcell{j} = X((j-1)*floor(n/k)+1:j*floor(n/k),:);
    end
    Xcell{k} = X((k-1)*floor(n/k)+1:n,:);
    
    for ii = 1:n_mc
        
        beta = normrnd(0,alpha/sqrt(p),[p,1]);
        
        %% search over lambdas
        lambda_fac = linspace(0.1,2,T);
        ni_lambda_grid = p/(alpha^2)*lambda_fac;
        
        MSE_usual_lambda = zeros(T,1); 
        for i=1:T
            ni_lambda = ni_lambda_grid(i);
            v = zeros(k,1);
            A = zeros(k,k);
            R = zeros(k,k);
            
            Qcell = cell(1,k);
            for j = 1:k
                Qcell{j} = (Xcell{j}'*Xcell{j}+ni_lambda*eye(p))^(-1)*(Xcell{j}'*Xcell{j});
            end
            
            for j = 1:k
                v(j) = beta'*Qcell{j}*beta;
            end
            
            for j = 1:k
                for u=1:k
                    A(j,u) = beta'*Qcell{j}*Qcell{u}*beta;
                end
            end
           
            for j = 1:k
                R(j,j) = trace((Xcell{j}'*Xcell{j}+ni_lambda*eye(p))^(-2)*Xcell{j}'*Xcell{j});
            end
            
        MSE_usual_lambda(i) = norm(beta)^2 - v'*(A+R)^(-1)*v;
         end
        
        MSE_usual = MSE_usual + MSE_usual_lambda;
    end
    
    MSE_mean = MSE_usual/n_mc;
    
    %%
    rng(2);
    savefigs=1;    closefigs=0; a = {'-','--','-.',':'};
    figure, hold on
    %mi = 0;
    %ma = max(max(MSE_mean), max(MSE_mean_n));
    %ylim([mi,ma]);
    
    h1 = plot(lambda_fac,MSE_mean,'LineWidth',3,'color',rand(1,3));
    set(h1,'LineStyle',a{1});
    pa=1; %your point goes here
    x=[pa,pa];
    y=get(gca,'Ylim');
    h2 = plot(x,y,'linewidth',3,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
    %h3 = plot(lambda_fac,MSE_mean_n,'LineWidth',3,'color',rand(1,3));
    %set(h3,'LineStyle',a{3});
    
    xlabel('lambda');
    ylabel('MSE');
    xlim([min(lambda_fac),max(lambda_fac)]);
    set(gca,'fontsize',20)
    legend([h1,h2],{'MSE','iso opt lambda'},'location','Best')
    %legend('Location','southwest','Theoretical','Numerical');
    %str = sprintf( 'alpha=%.2f-gamma=%.2f-m=%.d',alpha,gamma,m);
    %title(str);
    
    if savefigs==1
        filename = ...
            sprintf( './MSE-gen-cov-fn-lambda-design-%d-alpha=%.2f-gamma=%.2f-k=%.d-n_mc=%d.png',design,alpha,gamma,k,n_mc);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end
end
