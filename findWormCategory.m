function [leaveClusterLogInd, loneWormLogInd, inClusterLogInd,smallClusterLogInd] = ...
    findWormCategory(filename,inClusterNeighbourNum,minNeighbrDist,postExitDuration)

% function takes the path of the skeleton file and various worm
% classification variables and returns logical indices for leaveCluster and
% loneWorms according to the classificiation definitions.

%% INPUTS:
% filename: full path to the _skeletons.hdf5 file
% inClusterNeighbourNum: [1x1] double = 3. The number of close neighbors needed for worm to be considered 'in cluster'
% minNeighbrDist: [1x1] double = 2000. The minimum distance required for a worm to be considered 'lone worm'
% postExitDuration: [1x1] double = 5. The duration (in seconds) after a
% worm exits a cluster to be included in the leave cluster analysis. Input
% only needed if needing leaveClusterLogInd.

%% OUTPUTS:
% leaveClusterLogInd: logical index of leave cluster worms in the size of '/trajectories_data'
% loneWormLogInd: logical index of lone worms in the size of '/trajectories_data'
% inClusterLogInd: logical index of in cluster worms in the size of '/trajectories_data'
% smallClusterLogInd: logical index of small cluster worms in the size of '/trajectories_data'

%% load file
min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
neighbr_dist = h5read(filename,'/neighbr_distances');
frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
trajData = h5read(filename,'/trajectories_data');

%% classify worms
%% identify small cluster worms
smallClusterLogInd = (num_close_neighbrs==2 & neighbr_dist(:,3)>=minNeighbrDist)...
    |(num_close_neighbrs==3 & neighbr_dist(:,4)>=minNeighbrDist)...
    |(num_close_neighbrs==4 & neighbr_dist(:,5)>=minNeighbrDist);
%% identify in cluster worms
inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;
%% identify lone worms
loneWormLogInd = min_neighbr_dist>=minNeighbrDist;
%% identify leave cluster worms
if nargin <4
    leaveClusterLogInd = [];
else
    % find worm-frames where inCluster changes from true to false
    leaveClusterLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
    leaveClusterStart = find(leaveClusterLogInd);
    % loop through each exit event, retain frames for the specified duration after a worm exits cluster
    for exitCtr = 1:numel(leaveClusterStart)
        thisExitIdx = leaveClusterStart(exitCtr);
        wormIndex = trajData.worm_index_joined(thisExitIdx);
        % check for the number of frames that the same worm has beyond the point of cluster exit
        wormPathLength = nnz(trajData.worm_index_joined(thisExitIdx:end)==wormIndex);
        if wormPathLength>=postExitDuration*frameRate
            leaveClusterEnd = leaveClusterStart+postExitDuration*frameRate-1;
        else
            leaveClusterEnd = leaveClusterStart+wormPathLength-1;
        end % this also excludes movie segments with ending frames beyond highest frame number
        % go through each starting frame to generate logical index for leave cluster worms
        leaveClusterLogInd(leaveClusterStart(exitCtr):leaveClusterEnd(exitCtr))=true;
    end
    % exclude when worms move back into a cluster
    leaveClusterLogInd(inClusterLogInd)=false;
    % exclude worms that have become lone worm
    leaveClusterLogInd(loneWormLogInd)=false;
end