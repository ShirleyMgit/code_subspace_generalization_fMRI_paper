root = '/home/fs0/abaram/scratch/shirley/alon';
addpath(genpath(fullfile(root,'code')));
addpath(genpath('/home/fs0/abaram/scratch/MATLAB/spm12'))
% subjects: sub-01 to sub-28
subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];            
%     runGlm1(root,subjects{iSub},1,1,0,false)    
%     createSubjMask_noNanOrInfAfterGlm(root,subjects{iSub},'glm1')
%     createSl(root,subjects{iSub},'glm1',100)
%     makeSubspaceGenerMatrix(root,sub,glm,'L100')
%     normAndSmoothMat(root,subjects{iSub},'L100')
%     calcProjMatContrastsSubj_forPalm(root,subjects{iSub},'L100')
end
% 
% % Run parametric group stats in MAtlab - to get a quick result
% groupLevelTTest(root,subjects,'L100')
% 
% %Prepare and run permutation group stats in Palm - including MC correction
% stackSubjectsContrasts(root,subjects)

nPerm = '10000';
conName = 'hexOnHexMinusClustOnHex';
maskName =  'juelich_EC10';
runPalm(root,conName,nPerm,maskName)

conName = 'visual_sameStructSameStimMinusSameStructDiffStim';
maskName =  'harvardOxford_iLOC50';
runPalm(root,conName,nPerm,maskName)
