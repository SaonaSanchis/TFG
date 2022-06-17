function [] = ClusterMKR_DSM(F_data,idx_cluster,n_dims_regression,Features_in,KNN1_vec)
%% Regression on cluster modes


uClust = unique(idx_cluster);
n_clust = length(uClust);
CT = linspecer(n_clust,'divergent');


for i_c = 1:n_clust

    disp(['Cluster ',num2str(i_c)])
    disp(['N = ',num2str(sum(idx_cluster == i_c))])
    disp('----------------------')
    %% COMPUTE PCA
    temp_F_data = F_data(idx_cluster == i_c,:);
    F_data_centered = temp_F_data - repmat(mean(temp_F_data),size(temp_F_data,1),1);
    [V,~] = eig(cov(F_data_centered));
    V = (fliplr(V))';
    F_data_PCA{i_c} = (V*F_data_centered')';
    Features_in_aux = Features_in;
    for a =1:length(Features_in)
        Features_in_aux{a}.Feature = Features_in{a}.Feature(:,idx_cluster == i_c);
    end
    
    %% MKR
    F_star = Plot_Variability_descriptor(Features_in_aux,F_data_PCA{i_c},n_dims_regression ,40*ones(1,11));
    
    %% PLOT FLOWS
    flow_names = {'Aortic Valve','Mitral Valve'};
    figure(1)
    set(gcf, 'Position', get(0, 'Screensize'));
    colororder(CT)
    for count_feat = 1:2
        subplot(3,3,count_feat)
        hold on
        if i_c == n_clust
            grid on
        end
        if count_feat == 2
            plot(linspace(0,100,length(F_star{count_feat}(:,3))),F_star{count_feat}(:,3),'LineWidth',4)
            ylim([0 150])
        else
            plot(linspace(0,100,length(F_star{count_feat}(:,3))),F_star{count_feat}(:,3),'LineWidth',4,'DisplayName',['Cluster ',num2str(i_c)])
            legend show
            ylim([-5 350])
        end
        xlabel('Normalized Cycle')
        ylabel([flow_names{count_feat},newline,'Flow (cm/s)'])
        
    end

       
    %% PLOT REGIONAL STRAINS
    if i_c == 4
        disp('P')
    end
    septal_names = {'Septal Basal','Septal Mid','Septal Apical'};
    
    for count_feat = 1:3
        subplot(3,3,3+count_feat)
        hold on
        if i_c == n_clust
            grid on
        end
        plot(linspace(0,100,length(F_star{count_feat+ 3}(:,3))),F_star{count_feat+3}(:,3),'LineWidth',4)
        disp(Features_in{count_feat+4}.str)
        xlabel('Normalized Cycle')
        ylabel(['LV',' ',septal_names{count_feat},newline,'(%)'])
        ylim([-30 10])
    end
    
    
    disp('=========================')
    lateral_names = {'Lateral Apical','Lateral Mid','Lateral Basal'};
    
    for count_feat = 1:3
        subplot(3,3,6+count_feat)
        hold on
        if i_c == n_clust
            grid on
        end
        plot(linspace(0,100,length(F_star{count_feat+ 6}(:,3))),F_star{count_feat+ 6}(:,3),'LineWidth',4)
        xlabel('Normalized Cycle')
        ylabel(['LV',' ',lateral_names{count_feat},newline,'(%)'])
        
        ylim([-30 10])
        
    end
    
    
end
fig_name = ['/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/Figures/Clustering/KNN_',num2str(KNN1_vec),'/','CLUSTER_MEANS_',num2str(n_dims_regression),'_dims.png'];
mkdir (['/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/Figures/Clustering/KNN_',num2str(KNN1_vec),'/'])
saveas(gcf,fig_name)
close(gcf)

end
