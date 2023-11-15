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

plotSwarmChartFlag = false; 

% where to get the data from - only one should be true
peakFlag = false;
avgEffect0p01Mask = false;%true;
avgEffect0p05Mask = true;%false;


if peakFlag
    titleStr = 'peak';
    peakVoxIndsFSL = [31,58,16];
    peakVoxInds = peakVoxIndsFSL + 1; % Matlab indexing

else
    if avgEffect0p01Mask
        titleStr = 'avg p<0.01';
        mask = fullfile(root,'masks','mni','visual_sameStructSameStimMinusSameStructDiffStim_nPerm_10000_harvardOxford_iLOC50_vox_tstat_uncp_c1_maskT099.nii.gz');
    elseif avgEffect0p05Mask
        titleStr = 'avg p<0.05';        
        mask = fullfile(root,'masks','mni','visual_sameStructSameStimMinusSameStructDiffStim_nPerm_10000_harvardOxford_iLOC50_vox_tstat_uncp_c1_maskT095.nii.gz');
    end
    Vmask = spm_vol(mask);
    maskData = spm_read_vols(Vmask);
end
allElem = 1:16;
proHexElem = [1, 2, 5, 6];
proHexClElem = [9, 10, 13, 14];
proSameMapElem = [1, 6, 11, 16];
proSameStructElem = [2, 5, 12, 15];

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
proSameMap_allData = allData(proSameMapElem, :);
proSameStruct_allData = allData(proSameStructElem, :);

% proHex = mean(proHex_allData,1);
% proHexCl = mean(proHexCl_allData,1);
% for iSub=1:length(subjects)% I dont understand why Alon did it
%     toSubtract = mean(proHex(iSub));
%     proHex(iSub) = proHex(iSub) - toSubtract;
%     proHexCl(iSub) = proHexCl(iSub) - toSubtract;
% end

projMatAllSubj = reshape(allData,[4,4,length(subjects)]);

if plotSwarmChartFlag
    nSub=28;
    
    proHex = mean(proHex_allData, 1);
    proHexCl = mean(proHexCl_allData, 1);

    proSameMap = mean(proSameMap_allData, 1);
    proSameStruct= mean(proSameStruct_allData, 1);

    % Sample data
    plot_warm(1, proHex - proHexCl, 'HexOnHex - HexOnComm');
    plot_warm(2, proSameMap - proSameStruct, 'SameMap - SameStructure');


end

%%

hex11_data = squeeze(projMatAllSubj(1, 1, :));
cluster11_data = squeeze(projMatAllSubj(3, 3, :));

figure(20)
subplot(2,1,1)
plot(hex11_data, '.')
title('hex 11, mean ', mean(hex11_data))
subplot(2,1,2)
plot(cluster11_data, '.')
title('cluster 11, mean ', mean(cluster11_data))

[c, p ] = corrcoef(hex11_data, cluster11_data)



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

    figure(num_figure)
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