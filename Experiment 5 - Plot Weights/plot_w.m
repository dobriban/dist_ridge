%Plots of weights
cd('C:\Dropbox\Projects\Distributed Ridge\Experiments\Experiment 5 - Plot Weights')
addpath('C:\Dropbox\Projects\Distributed Ridge\Experiments\Code')

%%
alpha_arr = [0.5,1,4,8];
gamma = [0:0.1:10];
for k=1:5
    for i = 1:length(alpha_arr)
        alpha = alpha_arr(i);
        
        w  = opt_w(k,gamma,alpha);
        
        rng(2);
        figure, hold on
        savefigs=1;    closefigs=1; a = {'-','--','-.',':'};
        h1 = plot(gamma,w,'LineWidth',3,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        h2 = plot(gamma, 1/k*ones(1,length(gamma)),'LineWidth',3,'color',rand(1,3));
        set(h2,'LineStyle',a{2});
        
        xlabel('gamma');
        ylabel('weight');
        xlim([min(gamma),max(gamma)]);
        ylim([0,1]);
        set(gca,'fontsize',20)
        legend([h1,h2],{'opt weights','baseline weights'},'location','Best')
        
        
        str = sprintf( 'k=%d, alpha=%.2f',k, alpha);
        title(str);
        
        
        if savefigs==1
            filename = sprintf( './weight-plot-k=%d-alpha=%.2f.png',k,alpha);
            saveas(gcf, filename,'png');
            fprintf(['Saved Results to ' filename '\n']);
            if closefigs==1
                close(gcf)
            end
        end
    end
end

%% Surf+cont
alpha = linspace(0,5,30);
gamma = linspace(0,3,30);
[X,Y] = meshgrid(alpha, gamma);

for k=1:5
    %Z = are(X,Y);
    Z  = opt_w(k,Y,X);
    
    figure
    colormap(brighten(jet,0.6))
    surf(X,Y,Z)
    xlim([min(alpha),max(alpha)])
    ylim([min(gamma),max(gamma)])
    zlim([0,max(max(Z))])
    xlabel('\alpha');
    ylabel('\gamma');
    zlabel('W');
    %%
    savefigs=1;  closefigs=1;
    if savefigs==1
        filename = sprintf( './weight-surf-k=%d.png',k);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end
    
    %% Contour
    
    figure
    contour(X,Y,Z,'ShowText','on')
    xlim([min(alpha),max(alpha)])
    ylim([min(gamma),max(gamma)])
    xlabel('\alpha');
    ylabel('\gamma');
    colormap(brighten(jet,0.6))
    
    if savefigs==1
        filename = sprintf( './weight-cont-k=%d.png',k);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        if closefigs==1
            close(gcf)
        end
    end
end

