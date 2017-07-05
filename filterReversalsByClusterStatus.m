function [ clusterStatusReversalLogInd, interRevTimesFiltered, revDurationFiltered, incompleteInterRevFiltered ] = ...
    filterReversalsByClusterStatus(revStartInd, clusterStatusLogInd,...
    interRevTime, revDuration, incompleteInterRev)
% filters reversals by a given cluster status and adjusts interreversal
% times for when a cluster status changes between reversals
clusterStatusReversalLogInd = ismember(revStartInd,find(clusterStatusLogInd));
interRevTimesFiltered = interRevTime(clusterStatusReversalLogInd);
revDurationFiltered = revDuration(clusterStatusReversalLogInd);
% censor those inter-reversal times that arise from non-contiguous reversal
% sequences, eg when cluster status changes btw reversals
clusterStatusReversalInd = find(clusterStatusReversalLogInd);
clusterStatusInd = find(clusterStatusLogInd);
nonContRevs = find(diff(clusterStatusReversalInd)~=1);
for nCRctr = nonContRevs'
    % assign last time of same cluster status after each
    % non-contiguous reversal, and mark rev-intertime as censored
    interRevTimesFiltered(nCRctr) = ...
        clusterStatusInd(find(clusterStatusInd<revStartInd(clusterStatusReversalInd(nCRctr)+1),1,'last'))... % find last time of this cluster status that is smaller than the start of the next reversal
        -revStartInd(clusterStatusReversalInd(nCRctr)); % subtract the start of the current reversal
end
incompleteInterRevFiltered = incompleteInterRev(clusterStatusReversalLogInd);
incompleteInterRevFiltered(nonContRevs) = true;
end

