function runPalm(root,conName,nPerm,maskName)
% conName: name of contrast to run

% make sure palm is installed an in path.
palmDir = '/vols/Scratch/abaram/MATLAB/palm-alpha111/';
addpath(genpath(palmDir));

inputFile = fullfile(root,'subspaceGener','groupStats','stackedInputs',['con_' conName '.nii']);
outputDir = fullfile(root,'subspaceGener','groupStats','palm');
outputName = [conName '_nPerm_' nPerm '_' maskName];
maskFile = fullfile(root,'masks','mni',[maskName '.nii']);

palmStr = ['palm -i ' inputFile ' -T -n ' nPerm ' -o ' fullfile(outputDir,outputName) ' -ise -save1-p -m ' maskFile];
eval(palmStr)