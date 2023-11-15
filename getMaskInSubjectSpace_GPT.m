function getMaskInSubjectSpace_GPT(root, sub)

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)

maskName = 'EC_HexOnHexPeak_mask_100vox.nii';
maskDir = fullfile(root, 'masks', 'mni');
structural_dir = fullfile(root, 'struct', sub);
cst = cellstr(spm_select('FPList', structural_dir, '^iy.*'));
fnames = cellstr(spm_select('FPList', maskDir, '^EC_HexOnHex.*'));

% Resample the original mask to match subject space dimensions (106 x 106 x 72)
V = spm_vol(fnames{1});
Y = spm_read_vols(V);
VOX = [2 2 2]; % Voxel size of subject space (2mm isotropic)
bb = [NaN NaN NaN; NaN NaN NaN];
dims = [106 106 72];
matlabbatch{1}.spm.util.defs.comp{1}.def = cst;
matlabbatch{1}.spm.util.defs.out{1}.push.fnames = fnames;
matlabbatch{1}.spm.util.defs.out{1}.push.weight = '';
matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {maskDir};
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.bb = bb;
matlabbatch{1}.spm.util.defs.out{1}.push.fov.bbvox.vox = VOX;
matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];

% Clear previous matlabbatch and specify new normalization template
clear matlabbatch
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(cst);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(fullfile(maskDir, ['w' maskName]));
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = bb;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = VOX;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [sub '_'];

% Run the reslice and normalization jobs
spm_jobman('run', matlabbatch);
