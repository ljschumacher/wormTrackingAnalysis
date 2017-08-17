function [leaveClusterLogInd, loneWormLogInd] = findLeaveClusterWorms(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration)

% function takes the path of the skeleton file and various worm
% classification variables and returns logical indices for leaveCluster and
% loneWorms according to the classificiation definitions. 

%% INPUTS:
% filename: full path to the _skeletons.hdf5 file
% inClusterNeighbourNum: [1x1] double = 3. The number of close neighbors needed for worm to be considered 'in cluster'
% minNeighbrDist: [1x1] double = 2000. The minimum distance required for a worm to be considered 'lone worm'
% postExitDuration: [1x1] double = 5. The duration (in seconds) after a worm exits a cluster to be included in the leave cluster analysis

%% OUTPUTS:
% leaveClusterLogInd: logical index of leave cluster worms in the size of '/trajectories_data'
% loneWormLogInd: logical index of lone worms in the size of '/trajectories_data'

%% load file
min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
neighbr_dist = h5read(filename,'/neighbr_distances');
frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
trajData = h5read(filename,'/trajectories_data');

%% classify worms
% identify in cluster worms
inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
% find worm-frames where inCluster changes from true to false
leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end)); 
leaveClusterStart = find(leaveClusterLogInd);
% loop through each exit event, retain frames for the specified duration after a worm exits cluster
for exitCtr = 1:numel(leaveClusterStart)
    wormIndex = trajData.worm_index_joined(leaveClusterStart(exitCtr));
    % check for the number of frames that the same worm has beyond the point of cluster exit
    wormPathLength = numel(find(trajData.worm_index_joined==wormIndex));
    if wormPathLength>=postExitDuration*frameRate
        leaveClusterEnd = leaveClusterStart+postExitDuration*frameRate; 
    else
        leaveClusterEnd = leaveClusterStart+wormPathLength;
    end
end
% exclude movie segments with ending frames beyond highest frame number
leaveClusterEnd = leaveClusterEnd(leaveClusterEnd<=numel(leaveClusterLogInd)); 
% trim starting frame list accordingly, since the ending frame list may be shortened at the end of the movie
leaveClusterStart = leaveClusterStart(1:numel(leaveClusterEnd));
% go through each starting frame to generate logical index for leave cluster worms
for exitCtr = 1:numel(leaveClusterStart)
    leaveClusterLogInd(leaveClusterStart(exitCtr):leaveClusterEnd(exitCtr))=true;
end
% exclude when worms move back into a cluster
leaveClusterLogInd(inClusterLogInd)=false; 
% exclude worms that have become lone worm
loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
leaveClusterLogInd(loneWormLogInd)=false; 