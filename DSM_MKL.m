% THIS CODE TRAINS WITH HCMs AND THEIR RELATIVES AND PROJECTS THE REST OF
% THE COHORT.
clc;
clear all;
close all;

load('/Users/pablomarti/Documents/MATLAB/Research/Echo Reading and Analyis/DSM_data_processed.mat')

%% Definition of features
%%%%%%%%%%%%%
% ALL FEATS %
%%%%%%%%%%%%%
FEATURES{1,1} = cell2mat(cellfun(@(x) x.AorticValve_final, DATA, 'UniformOutput', 0)');
FEATURES{2,1} = cell2mat(cellfun(@(x) x.MitralValve_final, DATA, 'UniformOutput', 0)');
FEATURES{3,1} = cell2mat(cellfun(@(x) x.Events, DATA, 'UniformOutput', 0))';
for a = 1:6
    aux_feats_1 = cell2mat(cellfun(@(x) x.x4C_LV_RegionalStrain_aligned(:,a), DATA, 'UniformOutput', 0)');
    FEATURES{3+a,1} = aux_feats_1;
end
FEATURES{10,1} = cell2mat(cellfun(@(x) x.x4C_LV_GlobalStrain_aligned, DATA, 'UniformOutput', 0))';

%% MKL -- Unsupervised option ?- training
addpath /Users/pablomarti/Documents/MATLAB/Research/MKL
% Total number of subjects
N = numel(DATA);

% Vector with the feature kind for posterior assignment of the kernel type
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

Features{1}.Feature = FEATURES{1};
Features{1}.str = 'Aortic Valve';
Features{2}.Feature = FEATURES{2};
Features{2}.str = 'Mitral Valve';
Features{3}.Feature = FEATURES{3};
Features{3}.str = 'Events';
Features{4}.Feature = FEATURES{4};
Features{4}.str = 'LV Strain1';
Features{5}.Feature = FEATURES{5};
Features{5}.str = 'LV Strain2';
Features{6}.Feature = FEATURES{6};
Features{6}.str = 'LV Strain3';
Features{7}.Feature = FEATURES{7};
Features{7}.str = 'LV Strain4';
Features{8}.Feature = FEATURES{8};
Features{8}.str = 'LV Strain5';
Features{9}.Feature = FEATURES{9};
Features{9}.str = 'LV Strain6';
Features{10}.Feature = FEATURES{10};
Features{10}.str = 'LV Global Strain';

Features = Features';




save DSM_MKL

