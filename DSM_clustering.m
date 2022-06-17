clearvars;close all;clc;
%%
load('DSM_MKL.mat')
addpath '/Users/pablomarti/Documents/MATLAB/Research/linspecer'
addpath '/Users/pablomarti/Documents/MATLAB/Research/Tools'
addpath '/Users/pablomarti/Documents/MATLAB/Research/Regression_OutputSpaceDimensions'


%% Read clinical & outcome data
load('ClinTab_DSM.mat')
DSM_CLIN_TABLE  = DSM_CLIN_TABLE(:,[1:4,6:12,14,22,23,]);

ClinParams = str2double(table2array(DSM_CLIN_TABLE));
ClinParams = ClinParams(:,2:end);
Data_dictionary = DSM_CLIN_TABLE.Properties.VariableNames;
Data_dictionary = Data_dictionary(1,2:end)';

% Match with echo data
IDs_echo = cellfun(@(x) x.ID, DATA, 'UniformOutput', 0);
IDsClin_aux = cellstr(DSM_CLIN_TABLE.PatientId);
for a =1:numel(IDs_echo)
    IDs_echo{a} = IDs_echo{a}([1:end-4,end-1:end]);
end
[Lia,Locb] = ismember( IDs_echo,IDsClin_aux);
Locb(Locb==0) = [];
ClinParams = ClinParams(int64(Locb),:);

% Exclude clinical parameters with a unique value --- non-informative
to_exclude = false(size(ClinParams,2),1);
for i = 1:size(ClinParams,2)
    to_exclude(i) = numel(unique(ClinParams(:,i))) == 1;
end
Data_dictionary(to_exclude) = [];
ClinParams(:,to_exclude) = [];
ParamClin = ClinParams;
%%
for aux_knn = 1:length(KNN1_vec)
    
    F_data = F_data_struct{aux_knn};
    
    dimA = 1;
    dimB = 2;
    dimC = 3;
    limits(1) = min(F_data(:,dimA)) - 0.4*abs(min(F_data(:,dimA)));
    limits(2) = max(F_data(:,dimA)) + 0.4*abs(max(F_data(:,dimA)));
    limits(3) = min(F_data(:,dimB)) - 0.4*abs(min(F_data(:,dimB)));
    limits(4) = max(F_data(:,dimB)) + 0.4*abs(max(F_data(:,dimB)));
    limits(5) = min(F_data(:,dimC)) - 0.4*abs(min(F_data(:,dimC)));
    limits(6) = max(F_data(:,dimC)) + 0.4*abs(max(F_data(:,dimC)));
    
    
    for n_dims = 3:6
        X = F_data(:,1:n_dims);
        
        % Clustering
        % SILHOUETTE
        eva = evalclusters(X,'kmeans','silhouette','KList',[2:5]);
        n_clusters = eva.OptimalK;
        
        figure()
        bar(eva.CriterionValues)
        xlabel('Number of clusters','FontSize',15)
        ylabel('Silhouette value','FontSize',15)
        grid on
        
        idx_cluster = kmeans(X, n_clusters,'Replicate',2000);
        
        aux_clusters_vec = NaN(1,n_clusters);
        for aux_clusters = 1:n_clusters
            aux_clusters_vec(aux_clusters) = sum(idx_cluster==aux_clusters);
        end
        
        if min(aux_clusters_vec) ~= 1
            
            % Regression on clusters
            ClusterMKR_DSM(X,idx_cluster,n_dims,Features,KNN1_vec(aux_knn))
            
            CT = linspecer(n_clusters,'divergent');
            uClust = unique(idx_cluster);
            n_clusters = length(uClust);
            figure('units','normalized','outerposition',[0 0 1 1])
            CT = linspecer(n_clusters,'divergent');
            colororder(CT)
            for aux_subplot = 1:4
                subplot(2,2,aux_subplot)
                hold on
                grid on
                for a =1:n_clusters
                    xx = X(idx_cluster==a,1);
                    yy = X(idx_cluster==a,2);
                    zz = X(idx_cluster==a,3);
                    scatter3(xx , yy, zz,100,'filled','LineWidth',3,'DisplayName',['Cluster ',num2str(a),' (n=',num2str(sum(idx_cluster==a)),')'])
                end
                
                xlim([limits(1) limits(2)])
                ylim([limits(3) limits(4)])
                zlim([limits(5) limits(6)])
                if aux_subplot == 4
                    lgd =legend('location','best');
                    lgd.FontSize = 15;
                end
                xlabel('Dimension 1','FontSize',25)
                ylabel('Dimension 2','FontSize',25)
                zlabel('Dimension 3','FontSize',25)
                
                
                switch aux_subplot
                    case 1
                        view(3)
                    case 2
                        
                        view(0,0)
                    case 3
                        
                        view(0,90)
                    case 4
                        
                        view(90,0)
                end
            end
            sgtitle('Output Space','FontSize',35)
            
            fig_name = ['/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/Figures/Clustering/KNN_',num2str(KNN1_vec(aux_knn)),'/','CLUSTERING_NDIMS',num2str(n_dims),'.png'];
            saveas(gcf,fig_name)
            
            close(gcf)
            
            %%
            
            
            % Define row name as clinical variables
            table_stats.RowName = Data_dictionary;
            nClusters = numel(unique(idx_cluster));
            
            % Define column name as clusters
            table_stats.ColumnName = cell(nClusters,1);
            for i = 1:nClusters
                table_stats.ColumnName{1,i} = ['Cluster ' num2str(i),' (n=',num2str(sum(idx_cluster==i)),')'];
            end
            
            % Indicate the binary variables
            feat_type = zeros(1,size(ParamClin,2));
            for i = 1:size(ParamClin,2)
                cond1 = numel(unique(ParamClin(:,i)))==2;
                cond2 = any(isnan(ParamClin(:,i)));
                cond3= numel(unique(ParamClin(:,i)))==3;
                if cond1 || (cond2 && cond3)
                    feat_type(i) = 1;
                end
            end
            
            % Create variables
            countVar = zeros(sum(feat_type),nClusters);
            percentageVar = zeros(sum(feat_type),nClusters);
            meanVar = zeros(sum(~feat_type),nClusters);
            stdVar = zeros(sum(~feat_type),nClusters);
            counter1 = 1;
            counter2 = 1;
            nFeatures = size(ParamClin,2);
            
            data = cell(nFeatures,nClusters);
            
            Categorical_variables = false(1,size(ParamClin,2) );
            for i = 1:size(ParamClin,2)
                cat_aux(i) = numel(unique(ParamClin(~isnan(ParamClin(:,i)),i))) <= 5;
            end
            
            for ii = 1:nFeatures
                for jj = 1:nClusters
                    if ii == 101
                        disp('Hola')
                    end
                    if feat_type(ii) %% If feature is binary
                        countVar(counter1,jj) = nansum(ParamClin(idx_cluster==jj,ii));
                        percentageVar(counter1,jj) = 100 * countVar(counter1,jj) / sum(idx_cluster==jj);
                        data{ii,jj} = [num2str(countVar(counter1,jj)),' (',num2str(round(percentageVar(counter1,jj),2)),'%)'];
                    else            %% If feature is continous/categorical
                        if kstest(normalize(ParamClin(:,ii))) == 0
                            meanVar(counter2,jj) = mean(ParamClin(idx_cluster==jj,ii),'omitnan'); %% Omit NaN from calculations
                            stdVar(counter2,jj) = std(ParamClin(idx_cluster==jj,ii),'omitnan');   %% Omit NaN from calculations
                            data{ii,jj} = [num2str(round(meanVar(counter2,jj),2)),' ',char(177),' ',num2str(round(stdVar(counter2,jj),2))];
                        else
                            meanVar(counter2,jj) = median(ParamClin(idx_cluster==jj,ii),'omitnan'); %% Omit NaN from calculations
                            stdVar(counter2,jj) = std(ParamClin(idx_cluster==jj,ii),'omitnan');   %% Omit NaN from calculations
                            y_25 = quantile(ParamClin(idx_cluster==jj,ii),0.25);
                            y_75 = quantile(ParamClin(idx_cluster==jj,ii),0.75);
                            data{ii,jj} = [num2str(round(meanVar(counter2,jj),2)),' (',num2str(round(y_25,2)),' to ',num2str(round(y_75,2)),')'];
                        end
                    end
                end
                if feat_type(ii)
                    % Chi-square test to compare statistical significance among groups
                    % for a binary variable
                    [~,~,pval] = crosstab(idx_cluster,ParamClin(:,ii));
                    data{ii,nClusters+1} = pval;
                else
                    if cat_aux(ii) % Categorical vars with Chi square
                        [~,~,pval] = crosstab(idx_cluster,ParamClin(:,ii));
                        data{ii,nClusters+1} = pval;
                    else
                        % We first check wether the variable was sampled from a normal distribution or not
                        if kstest(normalize(ParamClin(:,ii))) == 1
                            data{ii,nClusters+1} = kruskalwallis(ParamClin(:,ii),idx_cluster,'off'); %% Display option = off
                        else
                            data{ii,nClusters+1} = anova1(ParamClin(:,ii),idx_cluster,'off');
                        end
                    end
                end
                
                if feat_type(ii)
                    counter1 = counter1+1;
                else
                    counter2 = counter2+1;
                end
            end
            table_stats.Data = data;
            data = [Data_dictionary,   data];
            
            xl_data_name = ['/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/Figures/Clustering/KNN_',num2str(KNN1_vec(aux_knn)),'/',num2str(n_dims),'_dims_clustering.csv'];
            
            %cell2csv(xl_data_name,data)
            
            table_stats.ColumnName{1,nClusters+1} = 'p-value';
            % Convert cell to a table and use first row as variable names
            T = cell2table(data(:,2:end),'VariableNames',table_stats.ColumnName(1,:),'RowNames',table_stats.RowName);
            
            % Write the table to a CSV file
            writetable(T,xl_data_name,'WriteRowNames',true)
            save (['Solution_KNN_',num2str(KNN1_vec(aux_knn)),'_',num2str(n_dims),'_dims'])
        end
    end
end