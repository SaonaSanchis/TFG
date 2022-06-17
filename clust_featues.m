for mtb = 1:49   
    b = subset(:,mtb);
    name = F_names(mtb);
    for n_dims = 1:5
                X = F_data(:,1:n_dims);
                age = DSM_CLIN_TABLE(:,3);
                age = table2array(age);
                gender = DSM_CLIN_TABLE(:,4);
                
                
                n_clusters = 5;
                idx_cluster = kmeans(X, n_clusters,'Replicate',2000);
                lb =['Cluster ' num2str(n_dims)];
                figure('name',title);
                hold on
            
                boxplot(b,idx_cluster);
                
                
                xlabel(lb)
                ylabel(name)
                
                hold off; 
            
            fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\F_data_boxplt_5clust_',num2str(mtb),'_',num2str(n_dims),'.png'];
                        saveas(gcf,fig_name)
                        
                        close(gcf)
     end
end

%%

% Categoical vaiables

for n_dims = 1:5
    X = F_data(:,1:n_dims);
    gender = DSM_CLIN_TABLE(:,4);
    gender= table2array(gender);

    idx_cluster = kmeans(X, n_clusters,'Replicate',2000);
    figure('name','Dimensionality reduced space provided by MKL')
    hold on
    n = [idx_cluster gender];

    bar(n)
    xlabel('Cluster')
    ylabel('Gender')
    legend('show');
    grid on;
 
    title('Output space');
    hold off; 

    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\F_data_barplot_',num2str(n_dims),'.png'];
            saveas(gcf,fig_name)
            
            close(gcf)
end

%%
for mtb = 1:94
    b = subset(:,mtb);
    name = F_names(mtb);
    for i = 1:1
        F_data_s = F_data_struct(i);
        F_data = cell2mat(F_data_s);
        figure('name','Dimensionality reduced space provided by MKL')
        hold on
        s=50;
        scatter3(F_data(:,1),F_data(:,2),F_data(:,3),s,b,'filled')
        colormap("winter")
        grid on;
        xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3'); 
        colorbar
        title(name,num2str(KNN1_vec(i)));
        hold off; 
    
    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\cmap_',num2str(mtb),'_',num2str(KNN1_vec(i)),'.png'];
                saveas(gcf,fig_name)
                
                close(gcf)
    end
end

%%

    b = subset20(:,605);
    name = "AGe";
   
    for i = 1:1
        F_data_s = F_data_struct(i);
        F_data = cell2mat(F_data_s);

        figure('name','Dimensionality reduced space provided by MKL')
        hold on
        s=50;
        scatter3(F_data(:,1),F_data(:,2),F_data(:,3),s,b,'filled')
        colormap("winter")
        grid on;
        xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3'); 
        colorbar
        title(name,num2str(KNN1_vec(i)));
        hold off; 
    
    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\cmap_AGE_',num2str(mtb),'_',num2str(KNN1_vec(i)),'.png'];
                saveas(gcf,fig_name)
                
                close(gcf)
    end


