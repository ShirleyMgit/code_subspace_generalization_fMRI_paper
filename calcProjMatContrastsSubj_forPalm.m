function calcProjMatContrastsSubj_forPalm(root,sub,slName)
% save single subject contrusts of Tim's PCA projections analysis
spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spm_path));
util_path = fullfile(root,'code','util');
addpath(genpath(util_path));

nFiles = 16;
smoothKernel = 6;
matFname = [slName '_projMat'];

projMat = zeros(91,109,91,nFiles); % dims of standard brain

subspaceGener_dir = fullfile(root,'subspaceGener',sub);
projMatFiles = cell(nFiles,1);
for f=1:nFiles
    projMatFiles{f} = [subspaceGener_dir,'/smth',num2str(smoothKernel) 'mni_' matFname num2str(f),'.nii'];
    V           = spm_vol(projMatFiles{f}); % header info      
    [projMat_currElement,~] = spm_read_vols(V,0);    
    projMat(:,:,:,f) = projMat_currElement;
end

subjContrastsDir = fullfile(subspaceGener_dir,'contrasts');
if ~exist(subjContrastsDir,'dir')
    mkdir(subjContrastsDir);
end

%% Calculate main Hex contrast - same labels as in the paper (e.g.
% HsCs denotes the element of the matrix corresponding to activity from the
% small hexagonal graph projected on eigenvectors calculated from the small
% (same stimuli-set) community-structure graph. 
% 1:4 : projections on eigenvectors from (3 averaged runs of) Hl.
% 5:8 : projections on eigenvectors from (3 averaged runs of) Hs.
% 9:12 : projections on eigenvectors from (3 averaged runs of) Cl.
% 13:16 : projections on eigenvectors from (3 averaged runs of) Cs.

% projection of Hex on hex: 1: HlHl 2: HsHl 5: HlHs 6: HsHs
proHex = squeeze(projMat(:,:,:,1) + projMat(:,:,:,2) + projMat(:,:,:,5) + projMat(:,:,:,6));
% projections of Hex on Cluster: 9: HlCl 10: HsCl 13: HlCs 14 HsCs
proHexCl = squeeze(projMat(:,:,:,9) + projMat(:,:,:,10) + projMat(:,:,:,13) + projMat(:,:,:,14) );
proHex_contrast = proHex - proHexCl; 
pathContrast = fullfile(subjContrastsDir,'con_hexOnHexMinusClustOnHex.nii');
save4Dnii(proHex_contrast,V,pathContrast)


%% Calculate stim set ID contrast
% HlHl, HsHs, ClCl, CsCs
projSameStructSameStim = squeeze(projMat(:,:,:,1) + projMat(:,:,:,6) + projMat(:,:,:,11) + projMat(:,:,:,16));
% HsHl, HlHs, CsCl, ClCs
projSameStructDiffStim = squeeze(projMat(:,:,:,2) + projMat(:,:,:,5) + projMat(:,:,:,12) + projMat(:,:,:,15));

projStimSet_contrast = projSameStructSameStim - projSameStructDiffStim;
pathContrast = fullfile(subjContrastsDir,'con_visual_sameStructSameStimMinusSameStructDiffStim.nii');
save4Dnii(projStimSet_contrast,V,pathContrast)






