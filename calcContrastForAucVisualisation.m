clear all

root = '/home/fs0/abaram/scratch/shirley/alon';
addpath(genpath(fullfile(root,'code')));
nVoxels = 100;
maskName = ['EC_HexOnHexPeak_mask_' num2str(nVoxels) 'vox'];
load(fullfile(root,'subspaceGener',['AUC_visualisation_' maskName '.mat']))

% subjects: sub-01 to sub-28
subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];
end
% nVoxels = size(betasToUse,3);

% which betas: (for later analysis)
nPiles= 10;
nMaps = 4;
nRuns = 4;

cumVar_allSubs = nan(nMaps,nMaps,nRuns,nVoxels,length(subjects)); % nMapToBeProjected (single run) x nMap to get EVs from (3 runs of) x ...
H_all = nan(16,length(subjects));
for iSub = 1: length(subjects)
    
    betasToUse = squeeze(betasAllSubj(:,:,:,iSub));
    betaMeanLOO = zeros(40,4,nVoxels);

    nr = 1:4;
    
    %% LOO averages
    for r = 1:4
        betaMeanLOO(:,r,:) = mean(betasToUse(:,nr~=r,:),2); % mean of the betas in 3 of the 4 runs
    end
    
    s1L = zeros(4,1);
    s2L = zeros(4,1);
    s3L = zeros(4,1);
    s4L = zeros(4,1);
    
    for r = 1:4
        
        % Do SVD on the data averaged from the 3 runs except current run.
        % U is the spatial components, S the scalings (eigenvalues), V the
        % temporal components.
        
        %hex 1
        hex1avLoo = squeeze(betaMeanLOO(1:10,r,:));
        hex1avLoo = hex1avLoo - repmat(mean(hex1avLoo),nPiles,1); % substract the mean over condition
        [U1,~,~] = svd(hex1avLoo'*hex1avLoo,'econ');%
        
        % hex 2
        hex2avLoo = squeeze(betaMeanLOO(11:20,r,:));
        hex2avLoo = hex2avLoo - repmat(mean(hex2avLoo),nPiles,1);
        [U2,~,~] = svd(hex2avLoo'*hex2avLoo,'econ');%
        
        % cluster 1:
        clus1avLoo = squeeze(betaMeanLOO(21:30,r,:));
        clus1avLoo = clus1avLoo - repmat(mean(clus1avLoo),nPiles,1);
        [U3,~,~] = svd(clus1avLoo'*clus1avLoo,'econ');%
        
        % cluster 2:
        clus2avLoo = squeeze(betaMeanLOO(31:40,r,:));
        clus2avLoo = clus2avLoo - repmat(mean(clus2avLoo),nPiles,1);
        [U4,~,~] = svd(clus2avLoo'*clus2avLoo,'econ');%
        
        % Calculate SVD on the data from the run itself.
        % Sr are the scalings (eigenvalues)
        
        %hex 1
        hex1run = squeeze(betasToUse(1:10,r,:));
        hex1run =  hex1run - repmat(mean(hex1run,1),nPiles,1);
        [~,S1r,~] = svd(hex1run'*hex1run,'econ');%
        
        %hex 2
        hex2run = squeeze(betasToUse(11:20,r,:));
        hex2run = hex2run - repmat(mean(hex2run,1),nPiles,1);
        [~,S2r,~] = svd(hex2run'*hex2run,'econ');
        
        %cluster 1:
        clus1run = squeeze(betasToUse(21:30,r,:));
        clus1run = clus1run - repmat(mean(clus1run),nPiles,1);
        [~,S3r,~] = svd(clus1run'*clus1run,'econ');
        
        %cluster 2:
        clus2run = squeeze(betasToUse(31:40,r,:));
        clus2run = clus2run - repmat(mean(clus2run,1),nPiles,1);
        [~,S4r,~] = svd(clus2run'*clus2run,'econ');
        
        %projections:
        % To hex 1:
        pr1 = zeros(4,nVoxels);
        % project to eigenvectors
        P11 = (hex1run*U1).^2;
        P21 = (hex2run*U1).^2;
        P31 = (clus1run*U1).^2;
        P41 = (clus2run*U1).^2;
        % sum variance explained
        pr1(1,:) = sum(P11)/trace(S1r);
        pr1(2,:) = sum(P21)/trace(S2r);
        pr1(3,:) = sum(P31)/trace(S3r);
        pr1(4,:) = sum(P41)/trace(S4r);
        % calculate area under the curve
        v1L = cumsum(pr1,2);
        cumVar_allSubs(:,1,r,:,iSub) = v1L;
        s1L(:,r)  = sum(v1L,2);
        
        % To hex 2:
        pr2 = zeros(4,nVoxels);
        % project to eigenvectors
        P12 = (hex1run*U2).^2;
        P22 = (hex2run*U2).^2;
        P32 = (clus1run*U2).^2;
        P42 = (clus2run*U2).^2;
        % sum variance explained
        pr2(1,:) = sum(P12)/trace(S1r);
        pr2(2,:) = sum(P22)/trace(S2r);
        pr2(3,:) = sum(P32)/trace(S3r);
        pr2(4,:) = sum(P42)/trace(S4r);
        % calculate area under the curve
        v2L = cumsum(pr2,2);
        cumVar_allSubs(:,2,r,:,iSub) = v2L;        
        s2L(:,r)  = sum(v2L,2);
        
        % To cluster 1:
        pr3 = zeros(4,nVoxels);
        % project to eigenvectors
        P13 = (hex1run*U3).^2;
        P23 = (hex2run*U3).^2;
        P33 = (clus1run*U3).^2;
        P43 = (clus2run*U3).^2;
        % sum variance explained
        pr3(1,:) = sum(P13)/trace(S1r);
        pr3(2,:) = sum(P23)/trace(S2r);
        pr3(3,:) = sum(P33)/trace(S3r);
        pr3(4,:) = sum(P43)/trace(S4r);
        % calculate area under the curve
        v3L = cumsum(pr3,2);
        cumVar_allSubs(:,3,r,:,iSub) = v3L;        
        s3L(:,r)  = sum(v3L,2);
        
        % To cluster 2:
        pr4 = zeros(4,nVoxels);
        % project to eigenvectors
        P14 = (hex1run*U4).^2;
        P24 = (hex2run*U4).^2;
        P34 = (clus1run*U4).^2;
        P44 = (clus2run*U4).^2;
        % sum variance explained, and normalise by the total variance in the
        % projected data
        pr4(1,:) = sum(P14)/trace(S1r);
        pr4(2,:) = sum(P24)/trace(S2r);
        pr4(3,:) = sum(P34)/trace(S3r);
        pr4(4,:) = sum(P44)/trace(S4r);
        % calculate area under the curve
        v4L = cumsum(pr4,2);
        cumVar_allSubs(:,4,r,:,iSub) = v4L;        
        s4L(:,r)  = sum(v4L,2);
        
    end    
    
    % get the average of the area under the curve across all runs
    ms1 = mean(s1L,2);
    ms2 = mean(s2L,2);
    ms3 = mean(s3L,2);
    ms4 = mean(s4L,2);
    
    H = zeros(1,16);
    H(1:4) = ms1; % projections on eigenvectors from Hex1 (actually from the mean betas of 3 runs of Hex 1)
    H(5:8) = ms2; % projections on eigenvectors from Hex2
    H(9:12) = ms3; % projections on eigenvectors from cluster1
    H(13:16) = ms4; % projections on eigenvectors from cluster2
    H_all(:,iSub) = H;
    clear H
end

cumVar_allSubs = squeeze(mean(cumVar_allSubs,3)); % average over runs)
AUC = squeeze(mean(cumVar_allSubs,3));
AUC_proHex = squeeze(AUC(1,1,:) + AUC(2,1,:) + AUC(1,2,:) + AUC(2,2,:));
AUC_proHexCl = squeeze(AUC(1,3,:) + AUC(2,3,:) + AUC(1,4,:) + AUC(2,4,:));
AUC_contrast  = AUC_proHex - AUC_proHexCl;
[h,p,~,stats] = ttest(AUC_contrast)

proHex_toPlot = squeeze(cumVar_allSubs(1,1,:,:) +  cumVar_allSubs(2,1,:,:) + cumVar_allSubs(1,2,:,:) + cumVar_allSubs(2,2,:,:));
proHexCl_toPlot = squeeze(cumVar_allSubs(1,3,:,:) +  cumVar_allSubs(2,3,:,:) + cumVar_allSubs(1,4,:,:) + cumVar_allSubs(2,4,:,:));
%% Plot first 10 PCs cumulative variance explained 

nSub=length(subjects);
figure
subplot(1,2,1)
hold on
shadedErrorBar(1:10,mean(proHex_toPlot(1:10,:),2),std(proHex_toPlot(1:10,:),false,2)/sqrt(length(subjects)),'lineprops','r')
shadedErrorBar(1:10,mean(proHexCl_toPlot(1:10,:),2),std(proHexCl_toPlot(1:10,:),false,2)/sqrt(length(subjects)),'lineprops','b')


subplot(1,2,2)
hold on
scatter(ones(nSub,1),AUC_proHex,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),AUC_proHexCl,3,'filled','b','MarkerFaceAlpha',0.5)
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], [mean(AUC_proHex),mean(AUC_proHex)],'r','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], [mean(AUC_proHexCl),mean(AUC_proHexCl)],'b','lineWidth',2)
errorbar(1+0.2,mean(AUC_proHex),std(AUC_proHex)/sqrt(nSub),'r')
errorbar(2-0.2,mean(AUC_proHexCl),std(AUC_proHexCl)/sqrt(nSub),'b')
for iSub=1:length(subjects)
    plot([1,2],[AUC_proHex(iSub),AUC_proHexCl(iSub)],'color',[0 0 0 0.5])
end
set(gca,'FontSize',14,'XTick',1:2,'XTickLabels',{'HexOnHex','HexOnComm'},'XTickLabelRotation', 45)


%% Plot just the difference 
figure 
hold on
scatter(zeros(nSub,1)+0.01*randn(nSub,1),AUC_contrast,10,'filled','r','MarkerFaceAlpha',0.5)
wBar = 0.2;
plot([-wBar, +wBar], [mean(AUC_contrast),mean(AUC_contrast)],'k','lineWidth',2)
errorbar(0,mean(AUC_contrast),std(AUC_contrast)/sqrt(nSub),'k')
xlim([-0.3,0.3])
plot([-0.3,0.3],[0,0],'k--')

%% This is another way to calculate the same - it is indeed the same result!

% Calculate main Hex contrast - same labels as in the paper (e.g.
% HsCs denotes the element of the matrix corresponding to activity from the
% small hexagonal graph projected on eigenvectors calculated from the small
% (same stimuli-set) community-structure graph. 
% 1:4 : projections on eigenvectors from (3 averaged runs of) Hl.
% 5:8 : projections on eigenvectors from (3 averaged runs of) Hs.
% 9:12 : projections on eigenvectors from (3 averaged runs of) Cl.
% 13:16 : projections on eigenvectors from (3 averaged runs of) Cs.

% projection of Hex on hex: 1: HlHl 2: HsHl 5: HlHs 6: HsHs
proHex = squeeze(H_all(1,:) + H_all(2,:) + H_all(5,:) + H_all(6,:));
% projections of Hex on Cluster: 9: HlCl 10: HsCl 13: HlCs 14 HsCs
proHexCl = squeeze(H_all(9,:) + H_all(10,:) + H_all(13,:) + H_all(14,:));
proHex_contrast = proHex - proHexCl; 
[h,p,~,stats] = ttest(proHex_contrast)

%% save to CSV file the DABEST likes
toCSV = cell(28*2,2);
for iRow = 1:28;
toCSV{iRow,1} = 'HexOnHex';
toCSV{iRow,2} = proHex(iRow);
toCSV{28 + iRow,1} = 'HexOnComm';
toCSV{28 + iRow,2} = proHexCl(iRow);
end

T = cell2table(toCSV,'VariableNames',{'Identifiers','Values'});
writetable(T,fullfile(root,'subspaceGener',['PairedTestVisualisation_' maskName '.csv']))


%% Plot paired difference - this is the figure we'll use
figure
dabest(fullfile(root,'subspaceGener',['PairedTestVisualisation_' maskName '.csv']),'Paired')

