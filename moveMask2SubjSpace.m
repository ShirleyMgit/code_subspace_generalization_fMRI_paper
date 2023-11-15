function moveMask2SubjSpace(root,sub,maskFile)
% maskFile would usually be 'EC_HexOnHexPeak_mask_100vox.nii'

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)

maskMni = fullfile(root,'masks','mni',maskFile);
structural_dir = fullfile(root,'struct',sub);
cst = cellstr(spm_select('FPList', structural_dir, '^iy.*')); % select files starting with "iy" - I think this is just the inverse deformation field from the registration

clear matlabbatch
    
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cst;
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {maskMni};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
    NaN NaN NaN];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [sub '_'];

spm_jobman('run', matlabbatch);

subjMaskFullPath = fullfile(root,'masks',sub,[sub '_' maskFile]);
movefile(fullfile(root,'masks','mni',[sub '_' maskFile]),subjMaskFullPath);


%% orient mask into functional space: that should be the NATIVE functional
cd (fullfile(root,'masks',sub))

Im1 = fullfile(root,'glms','glm1',sub,'beta_0001.nii'); % load any image in native functional space
Im2 = subjMaskFullPath;
fileName = subjMaskFullPath; % I think Im2 and fileName are the same in my case
safM = -1000;

matlabbatch{1}.spm.util.imcalc.input = {Im1 %make a mask from Im1
    Im2 %and apply to Im2
    };
matlabbatch{1}.spm.util.imcalc.output = fileName;
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression =  ['i2.*(i1>',num2str(safM),')'];
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch);