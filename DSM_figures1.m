clearvars;close all;clc;
%%
load('DSM_MKL.m')
addpath '/Users/pablomarti/Documents/MATLAB/Research/linspecer'
addpath '/Users/pablomarti/Documents/MATLAB/Research/Tools'
addpath '/Users/pablomarti/Documents/MATLAB/Research/Regression_OutputSpaceDimensions'

%% Read clinical & outcome data
load('ClinTab_DSM.mat')
DSM_CLIN_TABLE  = DSM_CLIN_TABLE(:,[1:4,6:14,22,23,]);

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
    

    %%
    for aux_feats = 1:length(Data_dictionary)
        
        Y = ClinParams(:,aux_feats);
        figure('Renderer', 'painters', 'Position', [10 10 1500 900])
        
        aux_Y = unique(Y(~isnan(Y)));
        
        if length(aux_Y)>10
            Y = discretize(Y,10);
            aux_Y = unique(Y);
            CT = linspecer(length(aux_Y),'sequential');
        else
            CT = linspecer(length(aux_Y),'qualitative');
        end
        
        
        for aux_fig = 1:4
            subplot(2,2,aux_fig)
            hold on
            grid on
            for a = 1:length(aux_Y)
                xx = F_data(Y == aux_Y(a),dimA);
                yy = F_data(Y == aux_Y(a),dimB);
                zz = F_data(Y == aux_Y(a),dimC);
                aux_color = CT(a,:);
                scatter3(xx , yy, zz,150,aux_color,'filled','MarkerEdgeColor','k','DisplayName',num2str(aux_Y(a)))
            end
            xlim([limits(1) limits(2)])
            ylim([limits(3) limits(4)])
            zlim([limits(5) limits(6)])
            
            xlabel('Dim 1')
            ylabel('Dim 2')
            zlabel('Dim 3')
            switch aux_fig
                case 1
                    view(3)
                case 2
                    view(0,0)
                    legend('Location','northeast')
                case 3
                    view(0,90)
                case 4
                    view(90,0)
            end
        end
        sgtitle(Data_dictionary{aux_feats})
        fig_name = ['/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/Figures/KNN_',num2str(KNN1_vec(aux_knn)),'/',Data_dictionary{aux_feats},'.png'];
        saveas(gcf,fig_name)
        close gcf
    end
end
