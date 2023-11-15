function normAndSmoothMat(root,sub,slName)

doSmooth=1;
doNorm=1;

smoothKernel = 6; % smoothing kernel

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)

structural_dir = fullfile(root,'struct',sub);
subspaceGener_dir = fullfile(root,'subspaceGener',sub);
cst = cellstr(spm_select('FPList', structural_dir, '^y.*')); % select files starting with "y" - I think this is just the deformation field from the registration
nFile = 16;
matFname = [slName '_projMat']; % file name of the 16 files of the 4x4 SVD projections matrix, without the number of the element

matFnames = cell(nFile,1);
for f=1:nFile
    disp( [subspaceGener_dir,'/',matFname,num2str(f),'.nii'])
    matFnames{f} = [subspaceGener_dir,'/',matFname,num2str(f),'.nii'];
end

clear matlabbatch
if doNorm
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cst;
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = matFnames;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [NaN NaN NaN
        NaN NaN NaN];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'mni_';
    
    spm_jobman('run', matlabbatch);
end


clear matlabbatch

if doSmooth
    
    %matlabbatch{1}.spm.spatial.smooth.data = [];
    matlabbatch{1}.spm.spatial.smooth.data =  cellstr(spm_select('FPList', subspaceGener_dir, '^mni.*'));%
    matlabbatch{1}.spm.spatial.smooth.fwhm   = [smoothKernel smoothKernel smoothKernel];
    matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{1}.spm.spatial.smooth.im     = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = ['smth',num2str(smoothKernel)];
    
    spm_jobman('run', matlabbatch)
end

