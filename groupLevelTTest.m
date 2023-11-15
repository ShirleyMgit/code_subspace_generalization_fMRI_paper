function groupLevelTTest(root,subjects,slName)
% get a tstat for group level using a t-test (for a quick parametric test -
% stats for the paper will be performed non-parametrically with PALM)

nFiles = 16;
smoothKernel = 6;
matFname = [slName '_projMat'];

niiDims = [91,109,91]; % dims of standard brain
projMat = zeros(niiDims(1),niiDims(2),niiDims(3),nFiles,length(subjects)); 
for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGener_dir = fullfile(root,'subspaceGener',sub);
    projMatFiles = cell(nFiles,1);
    for f=1:nFiles
        projMatFiles{f} = [subspaceGener_dir,'/smth',num2str(smoothKernel) 'mni_' matFname num2str(f),'.nii'];
        V           = spm_vol(projMatFiles{f}); % heatder info      
        [projMat_currElement,~] = spm_read_vols(V,0);    
        projMat(:,:,:,f,iSub) = projMat_currElement;
    end
end

% projection of Hex on hex: 1: HlHl 2: HsHl 5: HlHs 6: HsHs
projHex = squeeze(projMat(:,:,:,1,:) + projMat(:,:,:,2,:) + projMat(:,:,:,5,:) + projMat(:,:,:,6,:));
% projections of Hex on Cluster: 9: HlCl 10: HsCl 13: HlCs 14 HsCs
projHexCl = squeeze(projMat(:,:,:,9,:) + projMat(:,:,:,10,:) + projMat(:,:,:,13,:) + projMat(:,:,:,14,:) );
% projections of Cluster on Hex:
projClHex = squeeze(projMat(:,:,:,3,:) + projMat(:,:,:,4,:) + projMat(:,:,:,7,:) + projMat(:,:,:,8,:) );
% projections of Cluster on Cluster:
projClCl = squeeze(projMat(:,:,:,11,:) + projMat(:,:,:,12,:) + projMat(:,:,:,15,:) + projMat(:,:,:,16,:) );
% contrasts:
projHex_contrast = projHex - projHexCl; 
projHex_ClHex_contrast = projHex - projClHex; 
projClCl_ClHex_contrast = projClCl - projClHex; 
clear proHex proHexCl projClHex projClCl

%% Calculate stim set ID contrast
% HlHl, HsHs, ClCl, CsCs
projSameStructSameStim = squeeze(projMat(:,:,:,1,:) + projMat(:,:,:,6,:) + projMat(:,:,:,11,:) + projMat(:,:,:,16,:));
% HsHl, HlHs, CsCl, ClCs
projSameStructDiffStim = squeeze(projMat(:,:,:,2,:) + projMat(:,:,:,5,:) + projMat(:,:,:,12,:) + projMat(:,:,:,15,:));
projStimSet_contrast = projSameStructSameStim - projSameStructDiffStim;
clear projMat projSameStructSameStim projSameStructDiffStim

%% stats
% nVox = prod(niiDims);
% projHex_contrast_linVox = reshape(projHex_contrast,[nVox,length(subjects)]);
% tStats_linVox = zeros(nVox,1);
% for iVox = 1:nVox
%     [~,~,~,stats] = ttest(projHex_contrast_linVox(iVox,:),0,'tail','right');
%     tStats_linVox(iVox) = stats.tstat;
% end
% tStats = reshape(tStats_linVox,niiDims);
% save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_hexOnHexMinusHexOnCluster.nii'))
% 
% nVox = prod(niiDims);
% projHex_ClHex_contrast_linVox = reshape(projHex_ClHex_contrast,[nVox,length(subjects)]);
% tStats_linVox = zeros(nVox,1);
% for iVox = 1:nVox
%     [~,~,~,stats] = ttest(projHex_ClHex_contrast_linVox(iVox,:),0,'tail','right');
%     tStats_linVox(iVox) = stats.tstat;
% end
% tStats = reshape(tStats_linVox,niiDims);
% save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_hexOnHexMinusClusterOnHex.nii'))

nVox = prod(niiDims);
projClCl_ClHex_contrast_linVox = reshape(projClCl_ClHex_contrast,[nVox,length(subjects)]);
tStats_linVox = zeros(nVox,1);
for iVox = 1:nVox
    [~,~,~,stats] = ttest(projClCl_ClHex_contrast_linVox(iVox,:),0,'tail','right');
    tStats_linVox(iVox) = stats.tstat;
end
tStats = reshape(tStats_linVox,niiDims);
save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_ClOnClMinusClusterOnHex.nii'))

% projStimSet_contrast_linVox = reshape(projStimSet_contrast,[nVox,length(subjects)]);
% tStats_linVox = zeros(nVox,1);
% for iVox = 1:nVox
%     [~,~,~,stats] = ttest(projStimSet_contrast_linVox(iVox,:),0,'tail','right');
%     tStats_linVox(iVox) = stats.tstat;
% end
% tStats = reshape(tStats_linVox,niiDims);
% save4Dnii(tStats,V,fullfile(root,'subspaceGener','groupStats','tStat_visual_sameStructSameStimMinusSameStructDiffStim.nii'))

