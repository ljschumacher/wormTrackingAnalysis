function [ revStartClusterStatusLogInd, interRevTimesFiltered, revDurationFiltered, incompleteInterRevFiltered ] = ...
    filterReversalsByClusterStatus(revStartInd, clusterStatusLogInd,...
    interRevTime, revDuration, incompleteInterRev)
% filters reversals by a given cluster status and adjusts interreversal
% times for when a cluster status changes between reversals
revStartClusterStatusLogInd = ismember(revStartInd,find(clusterStatusLogInd)); % for each reversal, check if it's of the required cluster status
% for each reversal, find if the cluster status changes before the next
% reversal
revStartClusterStatusInd = find(revStartClusterStatusLogInd); 
notClusterStatusInd = find(~clusterStatusLogInd);
for revCtr = revStartClusterStatusInd'
    thisRevStart = revStartInd(revCtr);
    thisInterRevTime = interRevTime(revCtr);
    nextClusterStatusChange = notClusterStatusInd(find(notClusterStatusInd>thisRevStart,1,'first'));
    if nextClusterStatusChange<thisRevStart+thisInterRevTime
        % set end of reversal to when cluster status changes
        interRevTime(revCtr) = nextClusterStatusChange-thisRevStart;
        incompleteInterRev(revCtr) = true; % and mark as not fully tracked
    end
end
interRevTimesFiltered = interRevTime(revStartClusterStatusLogInd);
revDurationFiltered = revDuration(revStartClusterStatusLogInd);
incompleteInterRevFiltered = incompleteInterRev(revStartClusterStatusLogInd);
end

