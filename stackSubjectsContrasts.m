function stackSubjectsContrasts(root,subjects)

%  stack subjects contrasts. A single file of stacked contrasts is needed as
% inut to PALM. 

% get analyses names from example subject
exampleInputDir = fullfile(root,'subspaceGener','sub-01','contrasts');
cd(exampleInputDir)
tmpFname = dir('*con*'); 
fname_con = cell(size(tmpFname));
for iAnalysis = 1:length(fname_con)
    fname_con{iAnalysis} = tmpFname(iAnalysis).name(1:end-4); % get rid of .nii suffix
end

outputDir = fullfile(root,'subspaceGener','groupStats','palm');
mkdir(outputDir);

for iAnalysis = 1:length(fname_con)
    str = ['fslmerge -t ' fullfile(outputDir,fname_con{iAnalysis}) ' '];
    for iSub = 1:length(subjects)
        str = [str fullfile(root,'subspaceGener',subjects{iSub},'contrasts',[fname_con{iAnalysis} '.nii']) ' '];
    end
    system (str);
    gunzip(fullfile(outputDir,[fname_con{iAnalysis} '.nii.gz']));
    delete(fullfile(outputDir,[fname_con{iAnalysis} '.nii.gz']))
end

