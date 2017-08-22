function headAngleDiff = getHeadAngleDiff(xcoords, ycoords, marker, smoothing, frameRate)

% function returns a vector of head angle differences (in radian) between
% each frame of a given trajectory

%% INPUTS:
% xcoords: vector or matrix (if multiple nodes) containing x coordinates for a given trajectory.
% ycoords: vector or matrix (if multiple nodes) containing y coordinates for a given trajectory.
% marker: string, 'pharynx' or 'bodywall'.
% smoothing: logical, true or false.
% frameRate: [1x1] scalar, 3 or 9.
%% OUTPUT:
% headAngleDiff: vector of head angle differences (in radian) between each frame of the given trajectory

%% FUNCTION
% calculate head angles
[angleArray,meanAngles] = makeAngleArray(xcoords,ycoords);
headAngle = angleArray+meanAngles;
% take mean head angles for body wall marker
if strcmp(marker,'bodywall')
    headAngle = nanmean(headAngle(:,1:7),2); % take mean of the first 8 out of 49 nodes (head only)
end
% get angle difference
if smoothing
    smoothFactor = frameRate+1; % set smoothing to be over 1 second
    angleDiff = headAngle(2:end) - headAngle(1:end-1);
    totalSmoothedFrames = length(angleDiff)-smoothFactor;
    headAngleDiff = NaN(totalSmoothedFrames,1);
    for smoothFrameCtr = 1:totalSmoothedFrames
        headAngleDiff(smoothFrameCtr) = nanmean(angleDiff(smoothFrameCtr:smoothFrameCtr+smoothFactor));
    end
else
    headAngleDiff = headAngle(2:end) - headAngle(1:end-1);
end