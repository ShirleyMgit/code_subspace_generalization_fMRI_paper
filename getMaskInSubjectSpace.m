function getMaskInSubjectSpace(root,sub)


spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)

maskName = 'EC_HexOnHexPeak_mask_100vox.nii';
maskDir = fullfile(root,'masks','mni');
structural_dir = fullfile(root,'struct',sub);
cst = cellstr(spm_select('FPList', structural_dir, '^iy.*')); % select files starting with "iy" 
fnames = cellstr(spm_select('FPList', maskDir, '^EC_HexOnHex.*'));

clear matlabbatch

matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(cst);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fnames;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [nan nan nan; nan nan nan];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [sub '_'];

spm_jobman('run', matlabbatch);
