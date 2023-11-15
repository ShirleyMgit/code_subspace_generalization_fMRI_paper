function createSubjMask_noNanOrInfAfterGlm(root,sub,glm)
%%% create non-NaN non-inf mask in subject space

masksDir = fullfile(root,'masks',sub); 

% tmpFileName = fullfile(masksDir,[glm '_tmpToDel.nii']);
fileName = fullfile(masksDir,[glm '_noNanNoInf.nii']);
Im1 = fullfile(root,'glms', 'glm1',sub,'beta_0001.nii'); % load beta

matlabbatch{1}.spm.util.imcalc.input = {Im1
                                           };
matlabbatch{1}.spm.util.imcalc.output = fileName;
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression =  'i1<1000000';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch);

% clear matlabbatch
% 
% matlabbatch{1}.spm.util.imcalc.input = {tmpFileName
%                                                 };
% matlabbatch{1}.spm.util.imcalc.output = fileName;
% matlabbatch{1}.spm.util.imcalc.outdir = {''};
% matlabbatch{1}.spm.util.imcalc.expression =  'i1<1000000';
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
% 
% spm_jobman('run',matlabbatch);

% delete (tmpFileName)
