clear
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
%addpath(genpath(fullfile(root,'code')));
addpath(genpath("C:\Users\User\OneDrive\Desktop\spm12"))

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end

plotSwarmChartFlag = true; 

% where to get the data from - only one should be true
peakFlag = false;
avgEffect0p01Mask = false;
avgEffect0p05Mask = true;
EHR_julich = false;

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
    elseif EHR_julich
        titleStr = 'right EC masks';
        mask = fullfile(exp_root,'ROImasks','Entorhinal_R_juelich.nii');
    end
    Vmask = spm_vol(mask);
    maskData = spm_read_vols(Vmask);
end
allElem = 1:16;
proHexElem = [1, 2, 5, 6];
proHexClElem = [9, 10, 13, 14];
proClClElem = [11, 12, 15, 16];
proClHexElem = [3, 4, 7, 8];

allData = nan(length(allElem),length(subjects));

for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGenerDir = fullfile(root,'subspaceGener',sub);

    % Full matrix
    for iElem = 1:16
        fname{iElem} = fullfile(subspaceGenerDir,['smth6mni_L100_projMat' num2str(allElem(iElem)) '.nii']);
        V_all = spm_vol(fname{iElem});
        if peakFlag
            allData(iElem,iSub) = spm_sample_vol(V_all,peakVoxInds(1),peakVoxInds(2),peakVoxInds(3),0);
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
end

proHex_allData = allData(proHexElem,:);
proHexCl_allData = allData(proHexClElem,:);
proClCl_allData = allData(proClClElem, :);
proClHex_allData = allData(proClHexElem, :);


if plotSwarmChartFlag
    nSub=28;
    
    proHex = mean(proHex_allData, 1);
    proHexCl = mean(proHexCl_allData, 1);

    proClCl = mean(proClCl_allData, 1);
    proClHex= mean(proClHex_allData, 1);

    % Sample data
  %  plot_warm(1, proHex - proHexCl, 'HexOnHex - HexOnComm');
    plot_warm(2, proClCl - proHexCl, proClCl - proClHex, 'CommOnComm - diffStructure','HexOnComm','CommOnHex');


end

%%

function plot_warm(num_figure, data1, data2, ylabel_plot, data1_label, data2_label)
    % Calculate mean and standard deviation
    meanValue1 = mean(data1);
    stdValue1 = std(data1);

    meanValue2 = mean(data2);
    stdValue2 = std(data2);

    figure(num_figure)
    % Create a swarm plot
    swarmchart(ones(size(data1)),data1,'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
    
    hold on;

    swarmchart(3*ones(size(data2)),data2,'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
    
    hold on;
    % Plot mean value
    plot(1, meanValue1, 'r*', 'MarkerSize', 10);
    hold on;
    plot(3, meanValue2, 'r*', 'MarkerSize', 10);
    hold on;
    % Plot standard deviation range
    %xRange = [0.8, 1.2]; % x-axis range for std indicator
    %yRange = [meanValue - stdValue, meanValue + stdValue]; % y-axis range for std indicator
    errorbar(1,meanValue1,stdValue1/sqrt(length(data1)));
    errorbar(3,meanValue2,stdValue2/sqrt(length(data2)));
    xlim([-0.5,4.5])
    plot([-0.5,4.5],[0,0],'k--')
    hold off;
    
    % Adjust plot labels and legend
    xticks([1 ,3]);
    yticks();
    ylabel(ylabel_plot)
    xticklabels({data1_label, data2_label'})
end