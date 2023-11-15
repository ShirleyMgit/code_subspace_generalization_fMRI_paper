clear
close all
clc

root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
spm_path = 'C:\Program Files\MATLAB\R2023a\spm12';
addpath(genpath(spm_path));

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end

% where to get the data from - only one should be true
peakFlag = false;
avgEffect0p01Mask = true;
avgEffect0p05Mask = false;

if peakFlag
    titleStr = 'peak';
    peakVoxIndsFSL = [31,58,16];
    peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing

else
    if avgEffect0p01Mask
        titleStr = 'avg p<0.01';
%         mask = fullfile(root,'masks','mni','hexOnHexMinusClustOnHex_nPerm_10000_juelich_EC10_tfce_tstat_fwep_c1_mask_p01.nii');
        mask = fullfile(root,'masks','mni','hexOnHexMinusClustOnHex_nPerm_10000_juelich_EC10_vox_tstat_uncp_c1_mask_p01.nii');
    elseif avgEffect0p05Mask
        titleStr = 'avg p<0.05';        
%         mask = fullfile(root,'masks','mni','hexOnHexMinusClustOnHex_nPerm_10000_juelich_EC10_tfce_tstat_fwep_c1_mask_p05.nii');
        mask = fullfile(root,'masks','mni','hexOnHexMinusClustOnHex_nPerm_10000_juelich_EC10_vox_tstat_uncp_c1_mask_p05.nii');
    end
    Vmask = spm_vol(mask);
    maskData = spm_read_vols(Vmask);
end
allElem = 1:16;
proHexElem = [1,2,5,6];
proHexClElem = [9,10,13,14];
proClHexElem = [3, 4, 7, 8];
proClClElem = [11, 12, 15, 16];


proHex_allData = nan(length(proHexElem),length(subjects));
proHexCl_allData = nan(length(proHexClElem),length(subjects));
proClHex_allData = nan(length(proClHexElem),length(subjects));
proClCl_allData = nan(length(proClClElem),length(subjects));
allData = nan(length(allElem),length(subjects));

for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGenerDir = fullfile(root,'subspaceGener',sub);
    % Full matrix
    for iElem = 1:16
        fname{iElem} = fullfile(subspaceGenerDir,['smth6mni_L100_projMat' num2str(allElem(iElem)) '.nii']);
        V_all = spm_vol(fname{iElem});
        if peakFlag
            allData(iElem,iSub) = spm_sample_vol(V_all, peakVoxInds(1), peakVoxInds(2), peakVoxInds(3),0);
        else
            tmp = spm_read_vols(V_all);
            % mask by effect thresholded mask
            tmp = tmp(:);
            maskData = logical(maskData(:));            
            blob = tmp(maskData);
            % take the mean in th the mask
            allData(iElem,iSub) = mean(blob);             
        end
    end
    proHex_allData = allData(proHexElem,:);
    proHexCl_allData = allData(proHexClElem,:);
    proClHex_allData = allData(proClHexElem,:);
    proClCl_allData = allData(proClClElem,:);
end

projMatAllSubj = reshape(allData,[4,4,length(subjects)]);

% calculate contrasts within blob:
main_effect = mean(proHex_allData - proHexCl_allData);
hex_cl_hex = mean(proHex_allData - proClHex_allData);
cl_cl_hex = mean(proClCl_allData - proClHex_allData);

% stats on contrast:
[h_main, p_main, ~ , stat_main] = ttest(main_effect, 0, 'Tail', 'right');
[h_hex_cl_hex, p_hex_cl_hex, ~ , stat_hex_cl_hex] = ttest(hex_cl_hex, 0, 'Tail', 'right');
[h_cl_cl_hex, p_cl_cl_hex, ~ , stat_cl_cl_hex] = ttest(cl_cl_hex, 0);


