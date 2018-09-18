% The function is similar to findWormCatergory.m, except that it gives logical indices for 
% both cluster entering and exiting worms. Logical indices are also
% extended by given durations before AND after the event of cluster
% entry/exit.

function [enterClusterXLogInd,leaveClusterXLogInd,enterClusterStartLogInd,leaveClusterStartLogInd] = ...
    findClusterEntryExitX(filename,inClusterNeighbourNum,minNeighbrDist,preExitDuration,postExitDuration)

%% INPUTS:
% filename: string. Path to the skeletons.hdf5 file;
% inClusterNeighbourNum: [1x1] double, usually = 3. Number of close neighbors, used for defining in cluster worms. 
% minNeighbrDist: [1x1] double, usually = 2000. Minimum distance to the next neighbor, used in defining lone worms.
% preExitDuration: [1x1] double. Duration (in seconds) before a worm enters/exits a cluster to be included in the analysis,
% postExitDuration: [1x1] double. Duration (in seconds) after a worm enters/exits a cluster to be included in the analysis,

%% OUTPUTS:
% enterClusterXLogInd: logical index for entering cluster worms, with extended duration before/after point of entry.
% leaveClusterXLogInd: logical index for leaving cluster worms, with extended duration before/after point of exit.
% enterClusterStartLogInd: logical index for point of cluster entry.
% leaveClusterStartLogInd: logical index for point of cluster exit. 

%% FUNCTION:

% load data
min_neighbr_dist = h5read(filename,'/min_neighbr_dist');
num_close_neighbrs = h5read(filename,'/num_close_neighbrs');
trajData = h5read(filename,'/trajectories_data');
frameRate = double(h5readatt(filename,'/plate_worms','expected_fps'));
if ~isempty(strfind(filename,'r_X1_skeletons')) % check whether manually joined worm index exits
    trajDataWormIndex = trajData.worm_index_manual;
else
    trajDataWormIndex = trajData.worm_index_joined;
end

% identify in cluster worms
inClusterLogInd = num_close_neighbrs>=inClusterNeighbourNum;

% identify lone worms
loneWormLogInd = min_neighbr_dist>=minNeighbrDist;

%% identify cluster entry worms
% find worm-frames where inCluster changes from false to true
enterClusterXLogInd = vertcat(false,~inClusterLogInd(1:end-1)&inClusterLogInd(2:end));
enterClusterStartLogInd = enterClusterXLogInd;
enterClusterStart = find(enterClusterXLogInd);

% loop through each exit event, retain frames for the specified duration before/after a worm exits cluster
for entryCtr = 1:numel(enterClusterStart)
    thisEntryIdx = enterClusterStart(entryCtr);
    wormIndex = trajDataWormIndex(thisEntryIdx);
    
    % check for the number of frames that the same worm has before the point of cluster entry
    wormPathLengthBefore = nnz(trajDataWormIndex(1:thisEntryIdx)==wormIndex);
    if wormPathLengthBefore>=preExitDuration*frameRate
        enterClusterBeforeStart = enterClusterStart-preExitDuration*frameRate+1;
    else
        enterClusterBeforeStart = enterClusterStart-wormPathLengthBefore+1;
    end % this also excludes movie segments with starting frames beyond lowest frame number
    
    % check for the number of frames that the same worm has beyond the point of cluster exit
    wormPathLengthAfter = nnz(trajDataWormIndex(thisEntryIdx:end)==wormIndex);
    if wormPathLengthAfter>=postExitDuration*frameRate
        enterClusterAfterStart = enterClusterStart+postExitDuration*frameRate-1;
    else
        enterClusterAfterStart = enterClusterStart+wormPathLengthAfter-1;
    end % this also excludes movie segments with ending frames beyond highest frame number
    
    % go through each starting frame to generate logical index for leave cluster worms
    enterClusterXLogInd(enterClusterBeforeStart(entryCtr):enterClusterAfterStart(entryCtr))=true;
end

% exclude worms that have become lone worm
enterClusterXLogInd(loneWormLogInd)=false;

%% identify cluster exit worms
% find worm-frames where inCluster changes from true to false
leaveClusterXLogInd = vertcat(false,inClusterLogInd(1:end-1)&~inClusterLogInd(2:end));
leaveClusterStartLogInd = leaveClusterXLogInd;
leaveClusterStart = find(leaveClusterXLogInd);

% loop through each exit event, retain frames for the specified duration before/after a worm exits cluster
for exitCtr = 1:numel(leaveClusterStart)
    thisExitIdx = leaveClusterStart(exitCtr);
    wormIndex = trajDataWormIndex(thisExitIdx);
    
    % check for the number of frames that the same worm has before the point of cluster exit
    wormPathLengthBefore = nnz(trajDataWormIndex(1:thisExitIdx)==wormIndex);
    if wormPathLengthBefore>=preExitDuration*frameRate
        leaveClusterBeforeStart = leaveClusterStart-preExitDuration*frameRate+1;
    else
        leaveClusterBeforeStart = leaveClusterStart-wormPathLengthBefore+1;
    end % this also excludes movie segments with starting frames beyond lowest frame number
    
    % check for the number of frames that the same worm has beyond the point of cluster exit
    wormPathLengthAfter = nnz(trajDataWormIndex(thisExitIdx:end)==wormIndex);
    if wormPathLengthAfter>=postExitDuration*frameRate
        leaveClusterAfterStart = leaveClusterStart+postExitDuration*frameRate-1;
    else
        leaveClusterAfterStart = leaveClusterStart+wormPathLengthAfter-1;
    end % this also excludes movie segments with ending frames beyond highest frame number
    
    % go through each starting frame to generate logical index for leave cluster worms
    leaveClusterXLogInd(leaveClusterBeforeStart(exitCtr):leaveClusterAfterStart(exitCtr))=true;
end

% exclude worms that have become lone worm
leaveClusterXLogInd(loneWormLogInd)=false;

end