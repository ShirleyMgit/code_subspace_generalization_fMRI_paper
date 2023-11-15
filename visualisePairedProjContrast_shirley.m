clear
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
dataDir = fullfile(root,"beh","training");
%addpath(genpath(fullfile(root,'code')));
addpath(genpath("C:\Users\User\OneDrive\Desktop\spm12"))

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end

plotSwarmChartFlag = false; 

% where to get the data from - only one should be true
peakFlag = false;
avgEffect0p01Mask = false;
avgEffect0p05Mask = true;
alon_mask = false%true;
alon_peak = false%true;
EHR_julich = false;

if peakFlag
    if alon_peak
        titleStr = 'alon_peak';
        peakVoxIndsFSL = [32, 60,22];
        peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing
    else
        titleStr = 'peak';
        peakVoxIndsFSL = [31,58,16];
        peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing
    end

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
    elseif alon_mask
        titleStr = 'Alon effect';
        mask =  fullfile(root,'masks','mni', 'relationalStructure_smth5_rh_nPerm10000_clstrTh3p1_clusterm_mask.nii');
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

mean_matrixHex_suqres = mean(proHex_allData);
mean_matrixHexCl_suqres = mean(proHexCl_allData);
mean_matrixClCl_suqres = mean(proClCl_allData);
mean_matrixClHex_suqres = mean(proClHex_allData);

[h_hex_cl, p_hex_cl, ~, stats] = ttest(mean_matrixHex_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_cl_hex, p_cl_hex, ~, stats_clhex] = ttest(mean_matrixHex_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_hex, p_clcl_hex] = ttest(mean_matrixClCl_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_clcl_clhex, p_clcl_clhex] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_clhex_sym, p_clcl_clhex_sym] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres);


cl_contrast = mean_matrixClCl_suqres - mean_matrixClHex_suqres;
load(fullfile(dataDir,'pConCsCcor2.mat'))
load(fullfile(dataDir,'pConCsCncor2.mat'))
load(fullfile(dataDir,'pConCsCncor1.mat'))

cor_minus_notcor = pConCsCcor2 - pConCsCncor2;

[c, p] = corrcoef(cl_contrast, cor_minus_notcor);
[c1, p1] = corrcoef(cl_contrast, pConCsCncor1);
[c2, p2] = corrcoef(cl_contrast, pConCsCncor2);

if plotSwarmChartFlag
    nSub=28;
    
    proHex = mean(proHex_allData, 1);
    proHexCl = mean(proHexCl_allData, 1);

    proClCl = mean(proClCl_allData, 1);
    proClHex= mean(proClHex_allData, 1);

    % Sample data
    figure(100)
    plot_warm(1, proHex - proHexCl, 'HexOnHex - HexOnComm');
    title(['p_{val}=', num2str(p_hex_cl)])
    plot_warm(2, proHex - proClHex, 'HexOnHex - CommOnHex');
    title(['p_{val}=', num2str(p_cl_hex)])
    plot_warm(3, proClCl - proClHex, 'CommOnComm - CommOnHex');
    title(['p_{val}=', num2str(p_clcl_clhex)])
    plot_warm(4, proClCl - proHexCl, 'CommOnComm - HexOnComm');
    title(['p_{val}=',num2str(p_clcl_hex)])

end

%%

projMatAllSubj = reshape(allData,[4,4,length(subjects)]);

hex11_data = squeeze(projMatAllSubj(1, 1, :));
cluster11_data = squeeze(projMatAllSubj(3, 3, :));

figure(20)
subplot(2,1,1)
plot(hex11_data, '.')
title('hex 11, mean ', mean(hex11_data))
subplot(2,1,2)
plot(cluster11_data, '.')
title('cluster 11, mean ', mean(cluster11_data))

[c, p ] = corrcoef(hex11_data, cluster11_data);


projMatAllSubj2 = projMatAllSubj;
projMatAllSubj2(: ,4 ,:) = projMatAllSubj(:, 3, :);
projMatAllSubj2(: ,3 ,:) = projMatAllSubj(:, 4, :);
v = projMatAllSubj2(4, :, :);
projMatAllSubj2(4 ,: ,:) = projMatAllSubj2(3, :, :);
projMatAllSubj2(3 ,: ,:) = v;

figure
hold on
title([titleStr ', mean'])
projMatMean = squeeze(mean(projMatAllSubj2,3));
imagesc(projMatMean);
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

figure
hold on
title([titleStr ', t-value'])
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

figure(10)
subplot(1,2,1)
projMatMean = squeeze(mean(projMatAllSubj2,3));
imagesc(projMatMean);
title([titleStr ', mean'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

subplot(1,2,2)
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
title([titleStr ', t-value'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

function plot_warm(num_figure, data, ylabel_plot)
    % Calculate mean and standard deviation
    meanValue = mean(data);
    stdValue = std(data);

    subplot(2,2,num_figure)
    % Create a swarm plot
    swarmchart(ones(size(data)),data, 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
    
    hold on;
    
    % Plot mean value
    plot(1, meanValue, 'r*', 'MarkerSize', 10);
    
    % Plot standard deviation range
    xRange = [0.8, 1.2]; % x-axis range for std indicator
    yRange = [meanValue - stdValue, meanValue + stdValue]; % y-axis range for std indicator
    errorbar(1,meanValue,stdValue/sqrt(length(data)));
    xlim([-0.5,2.5])
    plot([-0.5,2.5],[0,0],'k--')
    hold off;
    
    % Adjust plot labels and legend
    xticks([]);
    yticks(0);
    ylabel(ylabel_plot)
end