function groupLevelPermStats(root)

% make sure palm is installed an in path.
palmDir = '/vols/Scratch/abaram/MATLAB/palm-alpha111/';
addpath(genpath(palmDir));

groupDir = fullfile(root,'subspaceGener','groupStats');

clusterThresh = 'None'; 
nPerm       = '10000';% number of permutations
maskName = 'EC';
groupLevelPermTests(rootData,'GLM3','xRun','lAccumbens20',clusterThresh,nPerm,maskName)

% mask to run the permutation tests in and perform multiple comparisons in.
% mask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces. This
% happened inside searchlightDefinitionSurfaceWrapper.m

inDir   = fullfile(groupDir,'stackedInputs');
cd(inDir)
fname = dir('*con*'); 

outDir = fullfile(groupDir,'perm',maskName);
mkdir(outDir);
maskFile = fullfile(rootData,'masks','mni',[maskName '.nii']);
for iAnalysis = 1:length(fname) % better to parallelise this
    outFile = fullfile(outDir,[fname(iAnalysis).name(1:end-4) '_nPerm' nPerm]);
    str = ['palm -i ' fname(iAnalysis).name ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile ' T'];      
    eval(str)
end
