clear all; close all; clc
%%
load DESMOPLAKYN_flows_strains_cropped.mat

%%
% FLOWS
exist_AV_PW = cell2mat(cellfun(@(x) isfield(x,'AorticValvePW_curves'), DATA, 'UniformOutput', 0));
exist_A_CW = cell2mat(cellfun(@(x) isfield(x,'AortaCW_curves') || isfield(x,'AorticValveCW_curves'), DATA, 'UniformOutput', 0)); % Only 1 has CW
exist_MV = cell2mat(cellfun(@(x) isfield(x,'MitralValvePW_curves'), DATA, 'UniformOutput', 0));
exist_TV = cell2mat(cellfun(@(x) isfield(x,'TricuspidValvePW_curves'), DATA, 'UniformOutput', 0));
exist_PV = cell2mat(cellfun(@(x) isfield(x,'PulmonaryValvePW_curves'), DATA, 'UniformOutput', 0));
% STRAINS
exist_LVStrain = cell2mat(cellfun(@(x) isfield(x,'x4CH_LV_GlobalStrain'), DATA, 'UniformOutput', 0));
% TISSUE VELOCITIES
exist_SeptalVel = cell2mat(cellfun(@(x) isfield(x,'SeptalMitralTDI_curves'), DATA, 'UniformOutput', 0));
exist_LateralVel = cell2mat(cellfun(@(x) isfield(x,'LateralMitralTDI_curves'), DATA, 'UniformOutput', 0));
exist_TricuspidVel = cell2mat(cellfun(@(x) isfield(x,'LateralTricuspidTDI_curves'), DATA, 'UniformOutput', 0));


ToKeep = exist_AV_PW & exist_A_CW & exist_MV & exist_TV & exist_LVStrain & exist_SeptalVel & exist_LateralVel & exist_TricuspidVel;
ToKeep_mat = [exist_AV_PW, exist_A_CW , exist_MV , exist_TV, exist_LVStrain , exist_SeptalVel , exist_LateralVel , exist_TricuspidVel,ToKeep];

ToKeep2 = exist_AV_PW & exist_MV  & exist_LVStrain;
ToKeep_mat2 = [exist_AV_PW, exist_MV, exist_LVStrain ,ToKeep2];

figure()
subplot(1,2,1)
imagesc(ToKeep_mat)
subplot(1,2,2)
imagesc(ToKeep_mat2)
%%
DATA = DATA(ToKeep2);
Subject = DATA;

%%
ECGs = cellfun(@(x) x.x4CH_LV_ECG, Subject, 'UniformOutput', 0);
StrainTimes = cellfun(@(x) x.x4CH_LV_StrainTimeAxis(2)-x.x4CH_LV_StrainTimeAxis(1), Subject, 'UniformOutput', 0);

figure()
hold on
for a = 1:length(ECGs)
    plot(StrainTimes{a},ECGs{a})
end


%%
Events_AVO = cell2mat(cellfun(@(x) x.AorticValvePW_valveOpening(1), Subject, 'UniformOutput', 0));
Events_AVC = cell2mat(cellfun(@(x) x.AorticValvePW_valveClosure(1), Subject, 'UniformOutput', 0));

Events_MVO = cell2mat(cellfun(@(x) x.MitralValvePW_valveOpening(1), Subject, 'UniformOutput', 0));
Events_MVC = cell2mat(cellfun(@(x) x.MitralValvePW_valveClosure(2), Subject, 'UniformOutput', 0));

Events_mat = [Events_AVO,Events_AVC,Events_MVO,Events_MVC];

figure()
hold on
for a = 1:4
    hist(Events_mat,25)
end

