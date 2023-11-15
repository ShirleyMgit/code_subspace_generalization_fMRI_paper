clear  

root = '/home/fs0/abaram/scratch/shirley/alon';
addpath(genpath(fullfile(root,'code')));
spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spm_path));
% subjects: sub-01 to sub-28

normAndRunGlm = false;

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end
nVoxels = 106;
% maskName = ['EC_HexOnHexPeak_mask_' num2str(nVoxels) 'vox'];
maskName = ['EC2_' num2str(nVoxels) 'vox'];
mask = fullfile(root,'masks','mni',[maskName '.nii']);

betasAllSubj = nan(40,4,nVoxels,28); % nPilesInRun x nRuns x nVox in SL x nSubj

for iSub=1:length(subjects)
    sub=subjects{iSub};
    if normAndRunGlm
        normRawData(root,sub)
        runGlm1(root,subjects{iSub},1,1,0,true)
    end
    %%
    load(fullfile(root,'glms','glm1_mni',subjects{iSub},'SPM.mat'))
    preproc_path = [root '/preproc/' subjects{iSub}];
    runs = {'run-01','run-02','run-03','run-04'};
    % loading images:
    inFiles = [];
    for iRun = 1:length(runs)
        run = runs{iRun};
        allfiles = cellstr(spm_select('FPList', [preproc_path,'/',run,'/'], '^mni_r_cleaned_smoothed_.*.nii$')); % pre warping!
        inFiles = [inFiles; allfiles];
    end
    VolIn=spm_vol(inFiles);
    X = nan(VolIn{1}.dim(1),VolIn{1}.dim(2),VolIn{1}.dim(3),length(VolIn));
    for iVol=1:length(VolIn)
        X(:,:,:,iVol) = spm_read_vols(VolIn{iVol});
    end
    Vmask = spm_vol({mask});
    maskData = spm_read_vols(Vmask{1});
    linVoxMask = find(maskData);
    [I,J,K] = ind2sub(size(maskData),linVoxMask);
    rawData = nan(length(linVoxMask),length(VolIn));
    for iVox=1:length(linVoxMask)
        rawData(iVox,:) = X(I(iVox),J(iVox),K(iVox),:);
    end
    [betaAll00,~,~,~,~,~]=rsa.spm.noiseNormalizeBeta(rawData',SPM);
    % which betas: (for later analysis)
    nPiles= 10;
    nMaps = 4;
    
    pilesBetasInds = ones(1,nPiles*nMaps);
    zeroBetas = 10 + 5 +2 +7; % 10: piles in map 5; 5: catch trials regs, one per block; 2: catch trial question and start experiment reg; 7: nuisance regressor (6 motion + CSF)
    betasToUse_singleRun =  [zeros(1,5) pilesBetasInds zeros(1,zeroBetas)]; % the first 5 0s: The regs modelling the random walks at the start of each map.
    betasToUse = [betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun];
    % get indeces of betas to use
    b_idx0 = 1:length(betasToUse);
    b_idx = b_idx0(betasToUse==1);
    
    % check there are voxels with no betas that are NaNs or Infs
    betasToUse = betaAll00(b_idx,:);
    isNanOrInf = isnan(mean(betasToUse,1))  + isinf(mean(betasToUse,1));
    
    if sum(isNanOrInf>0)
        disp(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
        warning(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
    end
    nVoxels = length(betasToUse(1,:));
    betasToUse = reshape(betasToUse,[40,4,nVoxels]); % nPiles*nMaps, nRuns, nVoxels
    betasAllSubj(:,:,:,iSub) = betasToUse;
end

save(fullfile(root,'subspaceGener',['AUC_visualisation_' maskName '.mat']),'betasAllSubj')
