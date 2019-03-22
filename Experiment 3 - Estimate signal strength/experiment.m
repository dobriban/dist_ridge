%Experiment with estimating signal strength
cd('C:\Dropbox\Projects\Distributed Ridge\Experiments\Experiment 3 - Estimate signal strength')
beep off
%% Experiment
rng(2);

n_grid = [50,500,5000,10000,50000];

rmse_sig = zeros(length(n_grid),1);
rmse_alpha = zeros(length(n_grid),1);

p = 50;
n_mc = 500;

for i=1:length(n_grid)
    i
    
    n = n_grid(i);
    e = zeros(n_mc,2);
    
    X = randn(n,p);
    [U,s,~] = svd(X);
    s_sq = diag(s).^2/p;
    if p<n
        s_sq = [s_sq; zeros(n-p,1)];
    end
    
    for j=1:n_mc
        j
        alpha = 1;
        beta = randn(p,1)*alpha/sqrt(p);
        
        sigma = 1;
        ep = randn(n,1)*sigma;
        y = X*beta+ ep;
        
        
        ty = U'*y;
        V = ty.^2;
        
        f = ones(n,2);
        f(:,1) = s_sq;
        e(j,:) = f\V;
    end
    
    err_sig = e(:,2)-sigma^2;
    err_alpha = e(:,1)-alpha^2;
    
    %% Plot err
    %rng(2);
    savefigs =1; %a = {'-','--','-.',':'};
    figure, hold on
    
    scatter(err_alpha,err_sig);
    xlabel('bias, alpha^2')
    ylabel('bias, sigma^2')
    set(gca,'fontsize',20)
    
    %legend(num2str(beta_arr'),'location','Best')
    
    if savefigs==1
        filename = sprintf( './est-sig-str-p=%d-n=%d-alpha=%.2f-n_mc=%d.png',p,n,alpha,n_mc);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
    
    %%
    rmse_sig(i) = sqrt(mean(err_sig.^2));
    rmse_alpha(i) = sqrt(mean(err_alpha.^2));
    
end

%%
rng(2);
savefigs =1; a = {'-','--','-.',':'};
figure, hold on
h1 = plot(n_grid,rmse_alpha,'linewidth',3,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(n_grid,rmse_sig,'linewidth',3,'color',rand(1,3));
set(h2,'LineStyle',a{2});
xlabel('n')
ylabel('RMSE')
set(gca,'fontsize',20)
xlim([min(n_grid), max(n_grid)]);


legend([h1,h2],{'RMSE \alpha','RMSE \sigma'},'location','Best')

if savefigs==1
    filename = sprintf( './rmse-n_mc=%d-p=%d-alpha=%.2f.png',n_mc,p,alpha); 
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    %close(gcf)
end




