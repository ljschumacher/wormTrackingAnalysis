function [ revFreq ] = countReversalFrequency( reversalLogInd, frameRate, ...
    signedSpeed, clusterStatusLogInd )
% counts the number of reversals and estimates a frequency
Nrev = nnz(reversalLogInd);
T = nnz(clusterStatusLogInd)/frameRate;
Trev = nnz(signedSpeed(clusterStatusLogInd)<0)/frameRate;
revFreq = Nrev./(T - Trev);

end

