
%Load data%

dat = readtable('dataPIE2.xlsx');
dat(1,:) = []  ;
prots = readtable('prots3.xlsx');
prots(1,:) = [] ;
DATA = innerjoin(dat,prots);
DATA(50,:) = [];

dat1 = table2cell(dat);
dats = dat1(:,2:403);
dats(50,:) =[];
dats = cell2mat(dats);

DSM_CLIN_TABLE = readtable('clinical.xlsx');
DSM_CLIN_TABLE  = DSM_CLIN_TABLE(1:105,[2 4 7:9]);
DSM_CLIN_TABLE(1,:) = [] ;
DSM_CLIN_TABLE.Var2{77,1} = 'PIE (5) 28CT';
DSM_CLIN_TABLE.Var2{78,1} = 'PIE (5) 29CT';
DSM_CLIN_TABLE.Properties.VariableNames([1]) = {'Var1'};
DATA_AND_CLINICAL = join(DATA,DSM_CLIN_TABLE);

subset1 = DATA_AND_CLINICAL(:,[4 6 9 13 15 20 21 27 32 36 51 52 81 90 108 139 155 184 193 204 216 225 135 138 347 381 406 408 409 413 419 420 421 426 427 428 429 431 433 449 453 454 458 461 467 489 542 552 560]);% 678 679 680 ]);
sub = table2cell(subset1);
subset2 = cell2mat(sub);
%% include preeclmpsia
%For all features

pre = xlsread('preclamsia.xlsx');
prep = xlsread('preclamsia.xlsx',2);
pcl = transpose(pre);
pclp = transpose(prep);
pclpm = [pcl,pclp];
pr = pclpm(:,1:675);
pcl_subset = pclpm(:,[4 6 9 13 15 20 21 27 32 36 51 52 81 90 108 139 155 184 193 204 216 225 135 138 347 381 406 408 409 413 419 420 421 426 427 428 429 431 433 449 453 454 458 461 467 489 542 552 560]);
subset3 = [subset2;pcl_subset];
subset = fillmissing(subset3, 'constant', 0);

%mtb1 = table2cell(DATA);
%mtb2 = mtb1(:,2:676);
%mtb = cell2mat(mtb2);
%subset1 = [mtb;pr];
%subset = fillmissing(subset1, 'constant', 0);
subset(isinf(subset)|isnan(subset)) = 0;
ss =[dats;pcl];

for i=1:49 % if we want to a the 3 clinical variables we change 49 to 51 
    FEATURES{i,1}=subset(:,i)';
end

F_names = subset1.Properties.VariableNames;


%%
subset = subset2;
subset(isinf(subset)|isnan(subset)) = 0;
subset = fillmissing(subset, 'constant', 0);
for i=1:52 % if we want to a the 3 clinical variables we change 49 to 51 
    FEATURES{i,1}=subset(:,i)';
end

F_names = subset1.Properties.VariableNames;

%% ALL variables minus clinical

d = DATA(:,2:676);
all_d = table2cell(d);
all_ds = cell2mat(all_d);
ALL = [all_ds;pr]; % all variables minus clinical
ALL(isinf(ALL)|isnan(ALL)) = 0;

%% MKL -- Unsupervised option ?- training
subset = ALL;
% Total number of subjects
N = numel(subset);

% Vector with the feature kind for posterior assignment of the kernael type
% and parameters
kind = ones(1,numel(FEATURES));

KNN1_vec = round(sqrt(N)/4:sqrt(N)/2:sqrt(N)*4);
for aux_knn = 1:numel(KNN1_vec)
    KNN1 = KNN1_vec(aux_knn);
    KNN2 = KNN1_vec(aux_knn);

    % Total number of features
    n_features = numel(FEATURES);

    for i = 1:n_features

        % Vector with the feature kind for posterior assignment of the kernel
        % type and parameters
        % kind = 1 --> Pattern
        % kind = 2 --> Continous variable
        % kind = 3 --> Binary variable
        % kind = 4 --> Categorical variable
        switch kind(i)
            case 1
                options.Kernel{i}.KernelType = 'exp_l2';
                options.Kernel{i}.Parameters = KNN1;
            case 2
                options.Kernel{i}.KernelType = 'exp_l2';
                options.Kernel{i}.Parameters = KNN1;
            case 3
                options.Kernel{i}.KernelType = 'binary';
                options.Kernel{i}.Parameters = N;
            case 4
                options.Kernel{i}.KernelType = 'ordinal';
                options.Kernel{i}.Parameters = N;
        end
    end

    options.AffinityNN = KNN2;
    
    % Unsupervised option
    [F_data_struct{aux_knn},betas{aux_knn},A{aux_knn}] = MKL(FEATURES,options);
end

for i = 1:50
   Features{i}.Feature = FEATURES{i};
  Features{i}.str = F_names(i);
end

Features = Features';

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

    %% Kernel correlation when increasing dimensions
    [F_data_self_corr,F_data_min_dims] = embedding_self_correlation(F_data);

    %%
    figure()
    hold on
    plot(F_data_self_corr,'LineWidth',5)
    % xline(F_data_min_dims,'--r',[num2str(F_data_min_dims),' Dims'],'LineWidth',.7,'HandleVisibility','off');
    yline(.9,'--r','90% Similarity','LineWidth',.7,'HandleVisibility','off');
    yline(1,'LineWidth',2,'HandleVisibility','off');
    scatter(1:length(F_data_self_corr),F_data_self_corr,150,'filled')

    ylim([-.05 1.2])
    xlim([0,50])
    grid on
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    ylabel('Kernel Alignment','FontSize',25)
    xlabel('Dimensions','FontSize',25)

    figure('units','normalized','outerposition',[0 0 1 1])
    for aux_subplot = 1:4
        subplot(2,2,aux_subplot)
        hold on
        grid on
        xx = F_data(:,1);
        yy = F_data(:,2);
        zz = F_data(:,3);
        scatter3(xx , yy, zz,100,'filled','LineWidth',3)
        xlim([limits(1) limits(2)])
        ylim([limits(3) limits(4)])
        zlim([limits(5) limits(6)])
        lgd =legend('location','best');
        lgd.FontSize = 15;
        xlabel('Dimension 1','FontSize',25)
        ylabel('Dimension 2','FontSize',25)
        zlabel('Dimension 3','FontSize',25)
        title(['KNN'])

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
    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\F_data_knn_',num2str(KNN1_vec(aux_knn)),'.png'];
    saveas(gcf,fig_name)

    close(gcf)
end

%% Dimensionality reduced space provided by MKL
%dats(50,:) =[]; % we use dats matrice when we don't include preclampsia
for i = 1:5
    F_data_s = F_data_struct(i);
    F_data = cell2mat(F_data_s);
    name = ["Hepatology","2","Pulmonary Disease","Preclampsia","Cardiac Disorder","Renal Condition"];
    figure('name','Dimensionality reduced space provided by MKL')
    hold on

    label = subset20(:,602);

    for iCluster = 1:6
        clustIdx =label==iCluster;
        plot3(F_data(clustIdx,1),F_data(clustIdx,2),F_data(clustIdx,3),'.','MarkerSize',30,...
            'DisplayName',sprintf(name(iCluster)));
    end

    legend('show');
    grid on;
    xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3');
    title('Output space');
    hold off;

    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\Condition_',num2str(KNN1_vec(i)),'.png'];
    saveas(gcf,fig_name)

    close(gcf)
end

%% Dimensionality reduced space provided by MKL
for i = 1:1
    F_data_s = F_data_struct(i);
    F_data = cell2mat(F_data_s);
    name2=["Healthy","Affected"];
    figure('name','Dimensionality reduced space provided by MKL')
    hold on
    label = subset20(:,603);
    for iCluster = 1:2
        clustIdx =label==iCluster;
        plot3(F_data(clustIdx,1),F_data(clustIdx,2),F_data(clustIdx,3),'.','MarkerSize',30,...
            'DisplayName',sprintf(name2(iCluster)));
    end
    legend('show');
    grid on;
    xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3');
    title('Output space');
    hold off;

    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\F_data_CC_',num2str(KNN1_vec(i)),'.png'];
    saveas(gcf,fig_name)

    close(gcf)
end
%% Dimensionality reduced space provided by MKL
for i = 1:1
    F_data_s = F_data_struct(i);
    F_data = cell2mat(F_data_s);
    name2=["Male","Feamle"];
    figure('name','Dimensionality reduced space provided by MKL')
    hold on
    label = subset20(:,604);
    for iCluster = 1:2
        clustIdx =label==iCluster;
        plot3(F_data(clustIdx,1),F_data(clustIdx,2),F_data(clustIdx,3),'.','MarkerSize',30,...
            'DisplayName',sprintf(name2(iCluster)));
    end
    legend('show');
    grid on;
    xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3');
    title('Output space');
    hold off;

    fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\F_data_Gender_',num2str(KNN1_vec(i)),'.png'];
    saveas(gcf,fig_name)

    close(gcf)
end


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


    for n_dims = 3:5
        X = F_data(:,1:n_dims);
        %age = DSM_CLIN_TABLE(:,3);
        %gender = DSM_CLIN_TABLE(:,4);
        % Clustering
        % SILHOUETTE
        %CC_cond = subset(:,1:2); % cojer las variables de sanos/affectados y tipo de condicion
        %X = [CC_cond,yy]; %juntar con la F_data
        eva = evalclusters(X,'kmeans','silhouette','KList',[2:8]);
        n_clusters = 5;

        figure()
        bar(2:8,eva.CriterionValues)
        xlabel('Number of clusters','FontSize',15)
        ylabel('Silhouette value','FontSize',15)

        grid on

        idx_cluster = kmeans(X, n_clusters,'Replicate',2000);

        %boxplot(idx_cluster,gender);
        %Case_control = tbl(:,2);
        %condition = tbl(:,3);
        %idx_cluster= condition;

        aux_clusters_vec = NaN(1,n_clusters);
        for aux_clusters = 1:n_clusters
            aux_clusters_vec(aux_clusters) = sum(idx_cluster==aux_clusters);
        end

        if min(aux_clusters_vec) ~= 1

            % Regression on clusters

            %ClusterMKR_DSM(X,idx_cluster,n_dims,Features,KNN1_vec(aux_knn))

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
            fig_name = ['C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\94_variables\kmenas' ...
                '_',num2str(KNN1_vec(aux_knn)),'_',num2str(n_dims),'.png'];
            saveas(gcf,fig_name)
            %close(gcf)

            %%

            %CC = ALL(:,1);
            %cc = (CC==1:2);
            %cond = ALL(:,2);
            %COND = (cond==1:6);
            %paramclin = [subset,cc];
            %parclin = [paramclin,COND];
            ParamClin=subset;
            Data_dictionary1=F_names;
            Data_dictionary =  Data_dictionary1.';
            % Define row name as clinical variables
            table_stats.RowName = Data_dictionary;
            nClusters = numel(unique(idx_cluster));
            %param_clin = input data
            %data_dictionary = col_names input data

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
                            data{ii,jj} = [num2str(round(meanVar(counter2,jj),2))];%,' ',char(177),' ',num2str(round(stdVar(counter2,jj),2))];
                        else
                            meanVar(counter2,jj) = median(ParamClin(idx_cluster==jj,ii),'omitnan'); %% Omit NaN from calculations
                            stdVar(counter2,jj) = std(ParamClin(idx_cluster==jj,ii),'omitnan');   %% Omit NaN from calculations
                            y_25 = quantile(ParamClin(idx_cluster==jj,ii),0.25);
                            y_75 = quantile(ParamClin(idx_cluster==jj,ii),0.75);
                            data{ii,jj} = [num2str(round(meanVar(counter2,jj),2))];%,' (',num2str(round(y_25,2)),' to ',num2str(round(y_75,2)),')'];
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



            mkdir ( strcat ('C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\KNN_94_var',num2str(KNN1_vec(aux_knn))));
            path = strcat ('C:\Users\User\Documents\BI3\TFG\Code\Code\Figures\KNN_94_var',num2str(KNN1_vec(aux_knn)));
            xl_data_name = [path,'/',num2str(n_dims),'_dims_clustering.csv'];

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

%% Correlation

corr=corrcoef(ALL);
corr=abs(corr);
[row,col,v] = find(corr < 0.7);


[ Avec, Ind ] = sort(corr(:),1,'descend');
max_values = Avec(1:1100);
[ ind_row, ind_col ] = ind2sub(size(corr),Ind(675:1100));

ALL(:,ind_row) = [];

ALL(isinf(ALL)|isnan(ALL)) = 0;
%%
for i=1:298
    FEATURES{i,1}=ALL(:,i)';
end

subset_298 = DATA(:,1:298);
F_names = subset_298.Properties.VariableNames;

for i = 1:298
    Features{i}.Feature = FEATURES{i};
    Features{i}.str = F_names(i);
end

Features = Features';

%%
corr_298=corrcoef(ALL);
corr_298=abs(corr_298);
ssort_298 = sort(corr_298);

[ Avec1, Ind1 ] = sort(corr_298(:),1,'descend');
max_values1 = Avec1(298:1100);
[ ind_row1, ind_col1 ] = ind2sub(size(corr_298),Ind1(298:1100));

%%
corr=corrcoef(ALL);
corr=abs(corr);

[ Avec, Ind ] = sort(corr(:),1,'ascend');
min_values = Avec(1:100);
[ ind_row, ind_col ] = ind2sub(size(corr),Ind(1:100));

ALL = ALL(:,ind_row);

ALL(isinf(ALL)|isnan(ALL)) = 0;
%%
for i=1:298
    FEATURES{i,1}=ALL(:,i)';
end

subset_298 = DATA(:,1:298);
F_names = subset_298.Properties.VariableNames;

for i = 1:298
    Features{i}.Feature = FEATURES{i};
    Features{i}.str = F_names(i);
end

Features = Features';