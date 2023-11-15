function normRawData(root,sub)

doSmooth=1;
doNorm=1;

smoothKernel = 6; % smoothing kernel

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)

runs = {'run-01','run-02','run-03','run-04'};
structural_dir = fullfile(root,'struct',sub);
for iRun=1:4
    preproc_dir = fullfile(root,'preproc',sub,runs{iRun});
    cst = cellstr(spm_select('FPList', structural_dir, '^y.*')); % select files starting with "y" - I think this is just the deformation field from the registration
    fnames = cellstr(spm_select('FPList', preproc_dir, '^r_cleaned_smoothed.*'));
    
    clear matlabbatch
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cst;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = fnames;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
        NaN NaN NaN];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'mni_';
       
    spm_jobman('run', matlabbatch);
    
end