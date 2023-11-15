function runGlm1(root,sub,specify,estimate,contrasts,mniFlag)

% mniFlag: true if running in MNI space. In this case currently only run
% the GLM in the EC HexOnHex peak voxel blob. 

pileDur = 1.4; % duration of pile (sequnce of images) presentation in seconds 
nPiles = 10; % number of piles per map
nRep = 4; % number of repeatitions for each pile
nCatch = 5; % number of catch trial in each block
nT = nRep*nPiles + nCatch; % number of trials of each block including catch trials. 
nBlocks = 5; % number of blocks in each run. one block per map. Only 4 of these were used - the ones of maps 1-4. 
nMaps = 5; % note that map 5 was actually a repitition of either map 1 or 2 alternately. these were excluded from further analyses. 
nCatchQ = nCatch + 1; % number of questions after catch trial piles in each block, including the extra question at the end of the block


preproc_path = [root '/preproc/' sub];

if mniFlag
    maskPath = fullfile(root,'masks','mni');
    mask = cellstr(spm_select('FPList', maskPath,'EC_HexOnHexPeak_mask_100vox.nii'));
else
    maskPath = fullfile(root,'masks',sub);
    mask = cellstr(spm_select('FPList', maskPath,'brain.nii'));
end
if mniFlag
    glm = 'glm1_mni'; % main glm for SVD analysis.
else
    glm = 'glm1'; % main glm for SVD analysis.
end
spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spm_path)
addpath(genpath(fullfile(root,'code')));

runs={'run-01','run-02','run-03','run-04'}; % run names (AKA sessions)
nRuns = length(runs);

if specify
    
    clear matlabbatch

    glmResultsDir = fullfile(root,'glms',glm,sub);
    if ~exist(glmResultsDir)
        mkdir(glmResultsDir)
    end
    for iRun=1:nRuns % sessions
        run = runs{iRun};
        load(fullfile(root,'beh',sub,[run,'.mat'])); % load behaviour files into a structure called RSA_maps
        
        mapsAllBlocks = RSA_maps(1).maps; % which maps used for each block (including "map number 5", though in practice map 5 was a repetition of either map 1 or 2 and was not used for analysis"
        
        % get onsets for (non-catch trials) piles and for catch trials
        [pilesNotCatchPresOnsets,pilesCatchPresOnsets] = getPilesOnsets(RSA_maps,nMaps,nBlocks,nT,nRep,nCatch,nPiles);  
                
        for iBlock = 1:nBlocks
            map = mapsAllBlocks(iBlock); % between 1-5
            
            switch map; % map in current block 
                
                % Regressors modelling the random walk at the start of each
                % block. These will not be used in later analyses. 
                case 1
                    onsets_pictureMap1   = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic; % onset picture == onset motor response except for first and last image
                    onsets_pictureMap1  = onsets_pictureMap1(onsets_pictureMap1 >0 );
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).name     = ['pictureMap', num2str(map)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(1).onset    = onsets_pictureMap1;
                case 2
                    onsets_pictureMap2   = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic; %  onset picture == onset motor response except for first and last image
                    onsets_pictureMap2  = onsets_pictureMap2(onsets_pictureMap2 >0 );
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(2).name     = ['pictureMap', num2str(map)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(2).onset    = onsets_pictureMap2;
                case 3
                    onsets_pictureMap3  = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic; %  % onset picture == onset motor response except for first and last image
                    onsets_pictureMap3  = onsets_pictureMap3(onsets_pictureMap3 >0 );
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(3).name     = ['pictureMap', num2str(map)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(3).onset    = onsets_pictureMap3;
                case 4
                    onsets_pictureMap4   = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic; % % onset picture == onset motor response except for first and last image
                    onsets_pictureMap4  = onsets_pictureMap4(onsets_pictureMap4 >0 );
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(4).name     = ['pictureMap', num2str(map)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(4).onset    = onsets_pictureMap4;
                case 5
                    if rem(iRun,2)==0 % map 1
                        onsets_pictureMap5_1   = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic;  % onset picture == onset motor response except for first and last image
                        onsets_pictureMap5_1  = onsets_pictureMap5_1(onsets_pictureMap5_1 >0 );
                        name5 = 'pictureMap5_1';
                        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(5).onset    = onsets_pictureMap5_1;
                    else % map 2
                        onsets_pictureMap5_2   = RSA_maps(iBlock).tmg.tpictureM(RSA_maps(iBlock).tmg.tpictureM~=0)-RSA_maps(1).tmg.tstart_pic;  % onset picture == onset motor response except for first and last image
                        onsets_pictureMap5_2  = onsets_pictureMap5_2(onsets_pictureMap5_2 >0 );
                        name5 = 'pictureMap5_2';
                        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(5).onset    = onsets_pictureMap5_2;
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(5).name     = name5;
                otherwise
                    error('wrong map number, Map repsentation switch!')
            end % end of Switch mapsAllBlocks(iBlock)            
            condition = iBlock; % conditions 1-5 are the presentations of picture in the random walk part of each block - not used for analysis
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).duration = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).tmod     = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).pmod     = struct('name',{}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).orth     = 0;                        
        end
                
        matlabbatch{1}.spm.stats.fmri_spec.dir            = {glmResultsDir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = 1.45;   % TR
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
        if mniFlag
            allfiles = cellstr(spm_select('FPList', [preproc_path,'/',run,'/'], '^mni_r_cleaned_smoothed_.*.nii$'));%not warp!            
        else
            allfiles = cellstr(spm_select('FPList', [preproc_path,'/',run,'/'], '^r_cleaned_smoothed_.*.nii$'));%not warp!
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).scans = allfiles; % select all scans - dummies are out during preprocess in fsl        
        
        % Calculate regressors of piles - these are the important
        % regressoors which wil lbe used in further analyses.         
        for iBlock = 1:nBlocks
            map = mapsAllBlocks(iBlock);            
            map_actual = RSA_maps(iBlock).pic.map; % difference betwen map_actual and map is only for "map 5" blocks: map actual has the map subjects actually saw (either map 1 or 2). These were excluded from analysis. 
            
            if map<5
                nameB = ['map',num2str(map),'Pile'];
            else
                nameB = ['secBmap',num2str(map_actual),'Pile'];
            end
            for n=1:10 % 10 piles in each map
                c = condition + (map-1)*10 + n; % condition==nBlock because so far there was a regressor for the random walk of each block. so c is between 6 and 
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).name     = [nameB,num2str(n)];
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).onset    = squeeze(pilesNotCatchPresOnsets(map,n,:));
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).duration = pileDur;
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).tmod     = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod     = struct('name',{}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).orth     = 0;
            end
        end
        condition = condition + nBlocks*10;
        
        for iBlock = 1:nBlocks
            map_actual = RSA_maps(iBlock).pic.map;
            map = mapsAllBlocks(iBlock);
            if map<5
                nameB = ['map',num2str(map),'catchPile'];
            else
                nameB = ['secBmap',num2str(map_actual),'catchPile'];
            end
            
            c = condition + map;
            
            % regressors for catch pile presentation (not the question)
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).name     = nameB;%[nameB,num2str(n)];
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).onset    = pilesCatchPresOnsets(map,:);
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).duration = pileDur;
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).tmod     = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).pmod     = struct('name',{}, 'param', {}, 'poly', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(c).orth     = 0;
            
        end
        
        % regressors for start trial and catch trials questions       
        onsetsStartPicture = zeros(1,nBlocks);
        catchQuestionOnsets = zeros(nBlocks*nCatchQ,1);
        catchResponseOnsets = zeros(nBlocks*nCatchQ,1); % onsets of catch trial responses
         for iBlock = 1:nBlocks
            onsetsStartPicture(iBlock) = RSA_maps(iBlock).tmg.tstart_pic - RSA_maps(1).tmg.tstart_pic;                                     
            onsetCatchQ = [RSA_maps(iBlock).tmg.tquestion(RSA_maps(iBlock).tmg.tquestion>0)-RSA_maps(1).tmg.tstart_pic; RSA_maps(iBlock).tmg.tquestionEndB-RSA_maps(1).tmg.tstart_pic];
            catchQuestionOnsets((iBlock-1)*nCatchQ+1:iBlock*nCatchQ)       = onsetCatchQ;%RSA_maps(blok).tmg.tquestion(RSA_maps(blok).tmg.tanswer>0)-RSA_maps(1).tmg.tstart_pic;           
            % timing of catch pile question response. this was not used as regressors. 
            onsetCatchR = [RSA_maps(iBlock).tmg.tanswer(RSA_maps(iBlock).tmg.tanswer>0)-RSA_maps(1).tmg.tstart_pic; RSA_maps(iBlock).tmg.tanswerEndB(RSA_maps(iBlock).tmg.tanswerEndB>0)-RSA_maps(1).tmg.tstart_pic];
            catchResponseOnsets((iBlock-1)*nCatchQ+1:iBlock*nCatchQ)  = onsetCatchR;            
        end 
        
        condition = condition + 6;   % condition + nBlock + 1     
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).name     = 'catch';
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).onset    = catchQuestionOnsets;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).duration = 0;%pileDur;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).pmod     = struct('name',{}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).orth     = 0;
        
        % Regressor for the "start experiment" message
        condition = 1+condition;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).name     = 'start_picture';
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).onset    = onsetsStartPicture;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).duration = 0;%pileDur;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).pmod     = struct('name',{}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).cond(condition).orth     = 0;
                
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi     = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).regress   = struct('name', {}, 'val', {});

        % add nuisance regressors - currently only th 6 motion paramaters
        % and the average csf signal. 
        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).multi_reg = cellstr(spm_select('FPList', [root,'/nuisanceRegs/' sub], [run '.mat$']));

        matlabbatch{1}.spm.stats.fmri_spec.sess(iRun).hpf       = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.00001;
%         matlabbatch{1}.spm.stats.fmri_spec.mask             = mask;
          
        % cellstr(spm_select('FPList', [spm_path,'\tpm'],
        % 'mask_ICV.nii$'));%is it for normalized brain?
        %matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
        
        
        disp(['------------ specifying 1st level ' sub ', ' run ,' done ------------'])
    end
    spm_jobman('run', matlabbatch)
end

disp(['------------ 1st level ' sub ' done ------------'])


%% estimation
if estimate
    
    clear matlabbatch
    
    currdir = glmResultsDir;
    
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[currdir,'/SPM.mat']};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('run', matlabbatch)
    disp(['------------estimate ' sub ' done ------------'])
    
end
%% contrasts
if contrasts
    
    clear matlabbatch
    currdir = glmResultsDir;
    
    matlabbatch{1}.spm.stats.con.spmmat = {[currdir,'\SPM.mat']};
    
    con = 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pictureMap';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [1 1 1 1 0];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pictureMapHexClu';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [1 1 -1 -1 0];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMapHex_Cluster';
    v = ones(1,10);
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v v -v -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';

    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMapSameStim';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v -v -v v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap1';
    v = ones(1,10);
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap2';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,15) v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap3';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,25) v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap4';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,35) v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap1_2';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap1_3';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v zeros(1,10) -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap2_3';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,15) v -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap1_4';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,5) v zeros(1,20) -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap2_4';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,15) v zeros(1,10) -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'pileMap3_4';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,25) v  -v];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
    
    con = con + 1;
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.name    = 'ctch';
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.weights = [zeros(1,45) 1];
    matlabbatch{1}.spm.stats.con.consess{con}.tcon.sessrep = 'repl';
 
    
    matlabbatch{1}.spm.stats.con.delete = 1;
    
    spm_jobman('run', matlabbatch)
    
    disp(['------------ contrasts' sub ' done ------------'])
end


disp(['------------ all done ' sub ' ------------'])