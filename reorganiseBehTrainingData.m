close all
clear
clc

dataDir = '/vols/Scratch/abaram/shirley/alon/beh/training';

% read tables
middleTable = readtable(fullfile(dataDir,'middleTable.csv')); %mM
pilesTable  = readtable(fullfile(dataDir,'pilesTable.csv')); %mP
closerTable = readtable(fullfile(dataDir,'closerTable.csv')); %mQ
navigTable = readtable(fullfile(dataDir,'navigTable.csv')); %mT

subjects = unique(middleTable.Var1); % list of all subject names in 'sub-XX' format
nSub = length(subjects);

% organise according to the 4 different tasks:
% is it in the middle, continue the pile, which is closer, navigation
D = struct();

D.middle.hex_day1_isCorrect = cell(length(subjects),1);
D.middle.hex_day1_ratioCorrect = zeros(length(subjects),1);
D.middle.hex_day2_isCorrect = cell(length(subjects),1);
D.middle.hex_day2_ratioCorrect = zeros(length(subjects),1);
D.middle.clus_day3_isCorrect = cell(length(subjects),1);
D.middle.clus_day3_ratioCorrect = zeros(length(subjects),1);
D.middle.clus_day4_isCorrect = cell(length(subjects),1);
D.middle.clus_day4_ratioCorrect = zeros(length(subjects),1);

D.closer.hex_day1_isCorrect = cell(length(subjects),1);
D.closer.hex_day1_ratioCorrect = zeros(length(subjects),1);
D.closer.hex_day2_isCorrect = cell(length(subjects),1);
D.closer.hex_day2_ratioCorrect = zeros(length(subjects),1);
D.closer.clus_day3_isCorrect = cell(length(subjects),1);
D.closer.clus_day3_ratioCorrect = zeros(length(subjects),1);
D.closer.clus_day4_isCorrect = cell(length(subjects),1);
D.closer.clus_day4_ratioCorrect = zeros(length(subjects),1);

D.piles.hex_day1_isCorrect = cell(length(subjects),1);
D.piles.hex_day1_ratioCorrect = zeros(length(subjects),1);
D.piles.hex_day2_isCorrect = cell(length(subjects),1);
D.piles.hex_day2_ratioCorrect = zeros(length(subjects),1);
D.piles.clus_day3_isCorrect = cell(length(subjects),1);
D.piles.clus_day3_ratioCorrect = zeros(length(subjects),1);
D.piles.clus_day4_isCorrect = cell(length(subjects),1);
D.piles.clus_day4_ratioCorrect = zeros(length(subjects),1);

D.navig.hex_day1_nSteps = zeros(length(subjects),7,3); % nSubs x nBlocks_day1 x nInitDists (ignore first block)
D.navig.hex_day2_nSteps = zeros(length(subjects),8,3); % nSubs x nBlocks_day2 x nInitDists
D.navig.clus_day3_nSteps = zeros(length(subjects),8,3); % nSubs x nBlocks_day3 x nInitDists
D.navig.clus_day4_nSteps = zeros(length(subjects),8,3); % nSubs x nBlocks_day4 x nInitDists




for iSub=1:nSub %s
    % middle
    rowsOfSub = strcmp(subjects{iSub},middleTable.Var1); %subVM
    blockNum = middleTable(rowsOfSub,:).Var2; %trialM: between 1 and 38, the relevant ones are 1:32, with 8 blocks per day
    isCorrect = middleTable(rowsOfSub,:).Var9;% isCorM
    % use logical indexing for block inds
    indBlocksHexDay1 = (blockNum > 0) & (blockNum <= 8);
    indBlocksHexDay2 = (blockNum > 8) & (blockNum <= 16);
    indBlocksClusDay3 = (blockNum > 16) & (blockNum <= 24);
    indBlocksClusDay4 = (blockNum > 24) & (blockNum <= 32);
    D.middle.hex_day1_isCorrect{iSub} = isCorrect(indBlocksHexDay1);
    D.middle.hex_day1_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay1))/length(isCorrect(indBlocksHexDay1));
    D.middle.hex_day2_isCorrect{iSub} = isCorrect(indBlocksHexDay2);
    D.middle.hex_day2_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay2))/length(isCorrect(indBlocksHexDay2));
    D.middle.clus_day3_isCorrect{iSub} = isCorrect(indBlocksClusDay3);
    D.middle.clus_day3_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay3))/length(isCorrect(indBlocksClusDay3));
    D.middle.clus_day4_isCorrect{iSub} = isCorrect(indBlocksClusDay4);
    D.middle.clus_day4_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay4))/length(isCorrect(indBlocksClusDay4));
    
    clear trialsOfSub blockNum isCorrect indBlocksHexDay1 indBlocksHexDay2 indBlocksClusDay3 indBlocksClusDay4
    
    % piles
    
    rowsOfSub = strcmp(subjects{iSub},pilesTable.Var1); %subVQ
    blockNum = pilesTable(rowsOfSub,:).Var2; %trialP
    isCorrect = pilesTable(rowsOfSub,:).Var5;% isCorP
    % using logical indexing for block inds
    indBlocksHexDay1 = (blockNum > 0) & (blockNum <= 8);
    indBlocksHexDay2 = (blockNum > 8) & (blockNum <= 16);
    indBlocksClusDay3 = (blockNum > 16) & (blockNum <= 24);
    indBlocksClusDay4 = (blockNum > 24) & (blockNum <= 32);
    D.piles.hex_day1_isCorrect{iSub} = isCorrect(indBlocksHexDay1);
    D.piles.hex_day1_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay1))/length(isCorrect(indBlocksHexDay1));
    D.piles.hex_day2_isCorrect{iSub} = isCorrect(indBlocksHexDay2);
    D.piles.hex_day2_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay2))/length(isCorrect(indBlocksHexDay2));
    D.piles.clus_day3_isCorrect{iSub} = isCorrect(indBlocksClusDay3);
    D.piles.clus_day3_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay3))/length(isCorrect(indBlocksClusDay3));
    D.piles.clus_day4_isCorrect{iSub} = isCorrect(indBlocksClusDay4);
    D.piles.clus_day4_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay4))/length(isCorrect(indBlocksClusDay4));
    
    clear trialsOfSub blockNum isCorrect indBlocksHexDay1 indBlocksHexDay2 indBlocksClusDay3 indBlocksClusDay4
    
    rowsOfSub = strcmp(subjects{iSub},closerTable.Var1); %subVQ
    blockNum = closerTable(rowsOfSub,:).Var2; %trialQ
    isCorrect = closerTable(rowsOfSub,:).Var8;% isCorQ
    indBlocksHexDay1 = (blockNum > 0) & (blockNum <= 8);
    indBlocksHexDay2 = (blockNum > 8) & (blockNum <= 16);
    indBlocksClusDay3 = (blockNum > 16) & (blockNum <= 24);
    indBlocksClusDay4 = (blockNum > 24) & (blockNum <= 32);    
    D.closer.hex_day1_isCorrect{iSub} = isCorrect(indBlocksHexDay1);
    D.closer.hex_day1_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay1))/length(isCorrect(indBlocksHexDay1));
    D.closer.hex_day2_isCorrect{iSub} = isCorrect(indBlocksHexDay2);
    D.closer.hex_day2_ratioCorrect(iSub) = sum(isCorrect(indBlocksHexDay2))/length(isCorrect(indBlocksHexDay2));
    D.closer.clus_day3_isCorrect{iSub} = isCorrect(indBlocksClusDay3);
    D.closer.clus_day3_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay3))/length(isCorrect(indBlocksClusDay3));
    D.closer.clus_day4_isCorrect{iSub} = isCorrect(indBlocksClusDay4);
    D.closer.clus_day4_ratioCorrect(iSub) = sum(isCorrect(indBlocksClusDay4))/length(isCorrect(indBlocksClusDay4));

    clear trialsOfSub blockNum isCorrect indBlocksHexDay1 indBlocksHexDay2 indBlocksClusDay3 indBlocksClusDay4

    % navigation
    % each row in the table is a single choice. there are multiple choices
    % within a trial (until you reach from start image to the target). There
    % are 3 trials in each block, with initial distances 2,3 and 4. 
    rowsOfSub = strcmp(subjects{iSub},navigTable.Var1); % subVT
    blockNum = navigTable(rowsOfSub,:).Var2; % trialT0
    choice = navigTable(rowsOfSub,:).Var7;
    initDist = navigTable(rowsOfSub,:).Var4; % ds0: the initial distance for the current trial - 2, 3, or 4
    % For navigation task we ignore trial 1 ont he first day because in that trial we gave full feedback.  
    indBlocksHexDay1 = (blockNum > 1) & (blockNum <= 8); 
    indBlocksHexDay2 = (blockNum > 8) & (blockNum <= 16);
    indBlocksClusDay3 = (blockNum > 16) & (blockNum <= 24);
    indBlocksClusDay4 = (blockNum > 24) & (blockNum <= 32); 

    blocks1 = blockNum(indBlocksHexDay1);
    blocks2 = blockNum(indBlocksHexDay2);
    blocks3 = blockNum(indBlocksClusDay3);
    blocks4 = blockNum(indBlocksClusDay4);
    
    % nSteps is the number of steps ignoring steps where participants
    % pressed enter to resample the possible answers. 
    

    [~, D.navig.hex_day1_nSteps(iSub,:,:)] = calc_nSteps_navig(blocks1, choice(indBlocksHexDay1), initDist(indBlocksHexDay1));
    [~, D.navig.hex_day2_nSteps(iSub,:,:)] = calc_nSteps_navig(blocks2, choice(indBlocksHexDay2), initDist(indBlocksHexDay2));
    [~, D.navig.clus_day3_nSteps(iSub,:,:)] = calc_nSteps_navig(blocks3, choice(indBlocksClusDay3), initDist(indBlocksClusDay3));
    [~, D.navig.clus_day4_nSteps(iSub,:,:)] = calc_nSteps_navig(blocks4, choice(indBlocksClusDay4), initDist(indBlocksClusDay4));

end


save(fullfile(dataDir,'allTrainingData.mat'),'D');
