clear all; close all; clc;
%% Read Aortic, Mitral and pulmonary flow data
disp('===================')
disp('Reading flow data...')
disp('===================')
% rawFolder = 'Data/Rocket_segmentations/';
rawFolder = '/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/out/';

% Retrieving the number of files as the elements on the list
% I obviate all files starting by '.'
listing = dir([rawFolder,'/*/*/*/*.csv']);
SampleNames = fullfile('', {listing.name});
SampleFolders = fullfile('', {listing.folder});

remove = cell2mat(cellfun(@(x) x(1)=='.', SampleNames, 'UniformOutput', 0));
SampleNames = SampleNames(~remove);
SampleFolders = SampleFolders(~remove);
num_samples = length(SampleNames);

DATA = cell(num_samples,1);

% print = input('Would you like to plot the curves? \n\n 1) Yes \n 2) No \n\n --> ');

% For every sample
for k = 1:num_samples
    aux_check = importfile_2022([SampleFolders{k},'/',SampleNames{k}]);
    path_splitted = strsplit(SampleNames{k},'_');
    
    patient_ID = path_splitted{1};
    DATA{k}.(genvarname('ID')) = patient_ID;
    
    if isfile([SampleFolders{k},'/',SampleNames{k}])
        
        aux_listing = dir([SampleFolders{k},'/',SampleNames{k}]);
        flowcurves = importfile_2022([SampleFolders{k},'/',SampleNames{k}]);
        flowcurves_headers = flowcurves(1,:);
        flowcurves = flowcurves(2:end,:);
        for aux_curves = 1:length(flowcurves(:,1))
            if strcmp(flowcurves(aux_curves,2),'curves') % Include velocities in column 4
                DATA{k}.(genvarname(append(flowcurves(aux_curves,1),'_',flowcurves(aux_curves,2)))) = [str2num(flowcurves(aux_curves,3))', str2num(flowcurves(aux_curves,4))'];
            else
                DATA{k}.(genvarname(append(flowcurves(aux_curves,1),'_',flowcurves(aux_curves,2)))) = [str2num(flowcurves(aux_curves,3))'];
            end
        end
        
    else
        disp('Warning: Non readable csv. Check for format errors')
        disp(['File is ',SampleNames{k}])
    end
    
end


%% Read strain data
% GLAUNES MIGHT NOT WORK IN POST SYSTOLIC MOTION, APICAL ROCKING, SEPTAL
% FLASH


% Retrieving the number of files as the elements on the list
% I obviate all files starting by '.'
rawFolder = '/Users/pablomarti/Documents/MATLAB/Research/AMISH_DESM/StrainExports_DSM/';
listing = dir(rawFolder);
SampleNames = fullfile('', {listing.name});
remove = cell2mat(cellfun(@(x) x(1)=='.', SampleNames, 'UniformOutput', 0));
SampleNames = SampleNames(~remove)';
remove = cell2mat(cellfun(@(x) endsWith(x,'.CSV'), SampleNames, 'UniformOutput', 0));
SampleNames = SampleNames(~remove)';
num_samples = length(SampleNames);

aux_count = 1;
aux_recheck = 1;
for a = 1:num_samples
    path_splitted = strsplit(SampleNames{a},'_');
    if a>2
        if ~strcmp(path_splitted{1},patient_ID)
            aux_count = aux_count + 1;
        end
    end
    img_view = path_splitted{7};
%     heart_chamber = path_splitted{3}(1:end-4);
    strain_type = [img_view,'_LV'];
    patient_ID = path_splitted{2};
    aux_ID = find(cell2mat(cellfun(@(x) strcmp(patient_ID,x.ID), DATA, 'UniformOutput', 0)));
    if isempty(aux_ID)
        disp(patient_ID)
    end
    if strcmp(img_view,'4CH') && ~isempty(aux_ID)
        tmp_data = importdata([rawFolder,SampleNames{a}(1:end-4),'.txt']);
        tmp_data = tmp_data.data;
        
        figure(1)
        plot(tmp_data(:,end-1))
        grid minor
        hold on
        plot(tmp_data(:,end).*0.05-abs(min(tmp_data(:,end-1)))-20)
        title([patient_ID,newline,strain_type])
        
        
        [strain_qrs,~] = ginput(2);
        strain_qrs = round(strain_qrs);
        global_strain = tmp_data(strain_qrs(1):strain_qrs(2),end-1);
        regional_strain = tmp_data(strain_qrs(1):strain_qrs(2),3:end-2);
        ecg_strain = tmp_data(strain_qrs(1):strain_qrs(2),end);
        time_strain =  tmp_data(strain_qrs(1):strain_qrs(2),1);
        
        
        if strain_qrs(1)>strain_qrs(2)
            strains_to_recheck{aux_recheck} = [patient_ID,'_',strain_type];
            figure(2)
            plot(global_strain)
            hold on
            plot(ecg_strain.*0.03-abs(min(tmp_data(:,end-1)))-20)
            grid minor
            title([patient_ID,newline,strain_type])
            fig_name = [patient_ID,'_',strain_type,'.png'];
            saveas(gcf,fig_name)
            pause(.2);
            close(gcf)
            aux_recheck = aux_recheck + 1;
        else
            DATA{aux_ID}.(genvarname([strain_type,'_','Strain time axis'])) = time_strain;
            DATA{aux_ID}.(genvarname([strain_type,'_','ECG'])) = ecg_strain;
            DATA{aux_ID}.(genvarname([strain_type,'_','Regional strain'])) = regional_strain;
            DATA{aux_ID}.(genvarname([strain_type,'_','Global strain'])) = global_strain;
        end
        clf
    end
        
    
    
end

save DESMOPLAKYN_flows_strains_cropped
% clearvars -except DATA
% save DESMOPLAKIN_cropped_strains

% save /Users/pablomarti/Documents/MATLAB/Research/ADOL/Echo_DATA.mat