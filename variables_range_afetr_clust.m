
F_data_s = F_data_struct(1);
F_data = cell2mat(F_data_s);
X = F_data(:,1:5);
n_clusters = 5;

idx_cluster = kmeans(X, n_clusters,'Replicate',2000);
label= idx_cluster;

for iCluster = 1:5
        clustIdx =label==iCluster;
        plot3(F_data(clustIdx,1),F_data(clustIdx,2),F_data(clustIdx,3),'.','MarkerSize',30,...
           'DisplayName',sprintf(name(iCluster)));
end
%%
for mtb = 1:94
    b = subset(:,mtb);
    name = F_names(mtb);

    lb =['Clusters '];
    %figure('name',title);
    hold on

    boxplot(b,idx_cluster);
    
    xlabel(lb)
    ylabel(name)
    
    hold off; 

    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\F_data_boxplt_5clust_',num2str(mtb),'.png'];
            saveas(gcf,fig_name)
            
    close(gcf)
 end