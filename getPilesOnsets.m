function [pilesOnsets, catchOnsets] = getPilesOnsets(B,nMaps,nBlocks,nT,nRep,nCatch,nPiles)
% B is the behavioural file, previously called RSAm (or RSA_maps)
% nT - number of trials including catch trials
% nMaps - 5

pilesOnsets = zeros(nMaps,nPiles,nRep); % Non-catch trials. Order: Hex 1, Hex 2, C1, C2, Hex B.  
catchOnsets = zeros(nMaps,nCatch);

maps = B(1).maps; % order of maps in current session (Alon: not sure why there are 5 vectors in ths which are the same, prob saved redundantly with a vec for each block which is why the 1st one is arbitrarily chosen)
for block=1:nBlocks
    iCatch = 1;
    onsets_pictureTrial0 = B(block).tmg.tpictureT(:,1)- B(1).tmg.tstart_pic;
    catchTrialsInds = B(block).pic.catchT; % trials numbers which are catch trials
    curMap = maps(block);
    pileId = B(block).pic.pileTask; % which pile was used in the current trial (including catch trials)
    counterRepsPerPile = zeros(1,nPiles); % vector that tracks the repetitions that each pile already had
    for iTrial=1:nT % including catch trials
        if sum(catchTrialsInds == iTrial-1) == 0 % if the PREVIOUS trial was NOT a catch trial. The "-1" is because the piles we are excluding are the ones that come after a catch trial question.
           counterRepsPerPile(pileId(iTrial)) = counterRepsPerPile(pileId(iTrial))+1;
           pilesOnsets(curMap,pileId(iTrial),counterRepsPerPile(pileId(iTrial))) = onsets_pictureTrial0(iTrial);
        else % If the previous trial was a catch trial
           catchOnsets(curMap,iCatch) = onsets_pictureTrial0(iTrial);
           iCatch = iCatch+1;
        end
    end
end
    