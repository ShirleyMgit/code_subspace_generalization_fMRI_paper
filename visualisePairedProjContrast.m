clear
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
%addpath(genpath(fullfile(root,'code')));
%spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
spm_path = 'C:\Program Files\MATLAB\R2023a\spm12';
%addpath(genpath('C:\Program Files\MATLAB\R2023a\spm12'))
addpath(genpath(spm_path));


subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end

plotSwarmChartFlag = true;%false; 

% where to get the data from - only one should be true
peakFlag = false;
avgEffect0p01Mask = false;
avgEffect0p05Mask = true;%false;
EHR_julich = false;%true;

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
proHexElem = [1,2,5,6];
proHexClElem = [9,10,13,14];

proHex_allData = nan(length(proHexElem),length(subjects));
proHexCl_allData = nan(length(proHexClElem),length(subjects));
allData = nan(length(allElem),length(subjects));




for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGenerDir = fullfile(root,'subspaceGener',sub);
% 
%     
%     for iElem = 1:length(proHexElem)
%         proHex_fname{iElem} = fullfile(subspaceGenerDir,['smth6mni_L100_projMat' num2str(proHexElem(iElem)) '.nii']);
%         proHexCl_fname{iElem} = fullfile(subspaceGenerDir,['smth6mni_L100_projMat' num2str(proHexClElem(iElem)) '.nii']);
%         V_hex = spm_vol(proHex_fname{iElem});
%         V_hexCl = spm_vol(proHexCl_fname{iElem});
% 
%         if peakFlag            
%             proHex_allData(iElem,iSub) = spm_sample_vol(V_hex,peakVoxInds(1),peakVoxInds(2),peakVoxInds(3),0);
%             proHexCl_allData(iElem,iSub) = spm_sample_vol(V_hexCl,peakVoxInds(1),peakVoxInds(2),peakVoxInds(3),0);
%         else 
%             tmp = spm_read_vols(V_hex);
%             % mask by effect thresholded mask
%             tmp = tmp .* maskData;
%             % take the mean in th the mask
%             proHex_allData(iElem,iSub) = mean(tmp(:)); 
%          
%             tmp = spm_read_vols(V_hexCl);
%             % mask by effect thresholded mask
%             tmp = tmp .* maskData;
%             % take the mean in th the mask
%             proHexCl_allData(iElem,iSub) = mean(tmp(:));             
%   
%         end
%     end
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
    proHex_allData = allData(proHexElem,:);
    proHexCl_allData = allData(proHexClElem,:);
end

proHex = mean(proHex_allData,1);
proHexCl = mean(proHexCl_allData,1);
for iSub=1:length(subjects)
    toSubtract = mean(proHex(iSub));
    proHex(iSub) = proHex(iSub) - toSubtract;
    proHexCl(iSub) = proHexCl(iSub) - toSubtract;
end

% %% save to CSV file the DABEST likes
% toCSV = cell(28*2,2);
% for iRow = 1:28;
%     toCSV{iRow,1} = 'HexOnComm';
%     toCSV{iRow,2} = proHexCl(iRow);
%     toCSV{iRow + 28,1} = 'HexOnHex';
%     toCSV{iRow + 28,2} = proHex(iRow);
% end
% 
% T = cell2table(toCSV,'VariableNames',{'Identifiers','Values'});
% writetable(T,fullfile(root,'subspaceGener','PairedTestVisualisation_FromSmth6Mni_EcPeak.csv'))
% 
% figure
% dabest(fullfile(root,'subspaceGener','PairedTestVisualisation_FromSmth6Mni_EcPeak.csv'),'Paired')



%%

if plotSwarmChartFlag
    nSub=28;
    contrast = proHex - proHexCl;
    proHex = mean(proHex_allData,1);
    proHexCl = mean(proHexCl_allData,1);
    figure
    % Sample data
    data = contrast;
    
    % Calculate mean and standard deviation
    meanValue = mean(data);
    stdValue = std(data);
    
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
    ylabel('HexOnHex - HexOnComm')

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

[c, p ] = corrcoef(hex11_data, cluster11_data)
figure
hold on
title([titleStr ', mean'])
projMatMean = squeeze(mean(projMatAllSubj,3));
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
[~,~,~,S] = ttest(projMatAllSubj,50,'dim',3);
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
projMatMean = squeeze(mean(projMatAllSubj,3));
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
[~,~,~,S] = ttest(projMatAllSubj,50,'dim',3);
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